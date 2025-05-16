#include <stdio.h>
#include <stdlib.h>
#include <string.h> // 用于 memset, memcpy 等内存操作函数
#include <GL/glut.h> // 包含 GLUT 库
#include <math.h>   // 用于 fmax, fmin

// 定义网格大小 N
// 流体模拟的有效区域是 N x N 个网格单元。
// Stam 的实现使用 (N+2) x (N+2) 的网格，最外层用于处理边界条件。
#define N 64

// 宏定义，用于将二维索引 (i, j) 映射到一维数组的索引。
// (i) 表示列索引，(j) 表示行索引。
// (N+2) 是每行的总单元格数，包括边界。
#define IX(i, j) ((i) + (N + 2) * (j))

// 遍历所有内部网格单元的宏，不包括边界层 (i, j 从 1 到 N)。
#define FOR_EACH_CELL for (int i = 1; i <= N; i++) { for (int j = 1; j <= N; j++) {
#define END_FOR }}

// 交换两个浮点型指针的宏。
// 在迭代求解或更新数组时，常用于交换当前值和前一时刻的值。
#define SWAP(x0, x) { float *tmp = x0; x0 = x; x = tmp; }

// 全局变量，用于存储流体状态的数组。
// u, v: 当前时刻的速度场分量 (u 为 x 方向，v 为 y 方向)
// u_prev, v_prev: 前一时刻的速度场分量，也用作临时缓冲区
// dens_r, dens_g, dens_b: 当前时刻的密度值 (红、绿、蓝分量)
// dens_prev_r, dens_prev_g, dens_prev_b: 前一时刻的密度值，也用作临时缓冲区
static float *u, *v, *u_prev, *v_prev;
static float *dens_r, *dens_g, *dens_b;
static float *dens_prev_r, *dens_prev_g, *dens_prev_b;
static int *is_obstacle; // 障碍物标记数组：0为流体，1为障碍物

// 模拟参数
static float dt = 0.1f; // 时间步长 (delta t)
static float visc = 0.000001f; // 粘度系数 (viscosity)，设为一个很小的值以减少扩散影响
static float diff = 0.000001f; // 密度扩散系数 (diffusion)，设为一个很小的值以减少扩散影响

// 窗口尺寸
static int win_width = 512;
static int win_height = 512;

// 鼠标拖拽信息
static int mouse_down[3]; // 鼠标按键状态
static int last_x, last_y; // 上一帧鼠标位置
static float current_draw_color[3] = {1.0f, 1.0f, 1.0f}; // 当前绘制颜色 (R, G, B), 默认白色

//--------------------------------------------------------------------------------------
// 函数声明
//--------------------------------------------------------------------------------------
static int allocate_data(void);
static void free_data(void);
static void add_source(int gridSize, float *x, float *s, float dt);
static void set_bnd(int gridSize, int b, float *x);
static void lin_solve(int gridSize, int b, float *x, float *x0, float a, float c);
static void diffuse(int gridSize, int b, float *x, float *x0, float diff, float dt);
static void advect(int gridSize, int b, float *d, float *d0, float *u, float *v, float dt);
static void project(int gridSize, float *u, float *v, float *p, float *div);
static void vel_step(int gridSize, float *u, float *v, float *u0, float *v0, float visc, float dt);
static void dens_step(int gridSize, float *x, float *x0, float *u, float *v, float diff, float dt);

// 障碍物相关函数
static void init_obstacles(void);
static void enforce_obstacle_boundaries(int gridSize, float *x, int b);

// GLUT 相关函数
static void display(void);
static void idle(void);
static void reshape(int w, int h);
static void keyboard(unsigned char key, int x, int y);
static void mouse(int button, int state, int x, int y);
static void motion(int x, int y);

// 模拟步进函数
void simulate_step();
void init_simulation();


//--------------------------------------------------------------------------------------
// 函数实现
//--------------------------------------------------------------------------------------

/**
 * @brief 分配流体模拟所需的数据空间。
 * @return 1 成功，0 失败。
 */
static int allocate_data(void) {
    int size = (N + 2) * (N + 2); // (N+2)x(N+2) 的网格总大小
    u = (float *)malloc(size * sizeof(float));
    v = (float *)malloc(size * sizeof(float));
    u_prev = (float *)malloc(size * sizeof(float));
    v_prev = (float *)malloc(size * sizeof(float));
    
    // 为 RGB 密度通道分配内存
    dens_r = (float *)malloc(size * sizeof(float));
    dens_g = (float *)malloc(size * sizeof(float));
    dens_b = (float *)malloc(size * sizeof(float));
    dens_prev_r = (float *)malloc(size * sizeof(float));
    dens_prev_g = (float *)malloc(size * sizeof(float));
    dens_prev_b = (float *)malloc(size * sizeof(float));

    // 为障碍物标记数组分配内存
    is_obstacle = (int *)malloc(size * sizeof(int));

    if (!u || !v || !u_prev || !v_prev || !dens_r || !dens_g || !dens_b ||
        !dens_prev_r || !dens_prev_g || !dens_prev_b || !is_obstacle) {
        fprintf(stderr, "无法分配数据内存\n");
        return 0;
    }

    // 初始化所有数据为0
    memset(u, 0, size * sizeof(float));
    memset(v, 0, size * sizeof(float));
    memset(u_prev, 0, size * sizeof(float));
    memset(v_prev, 0, size * sizeof(float));
    
    memset(dens_r, 0, size * sizeof(float));
    memset(dens_g, 0, size * sizeof(float));
    memset(dens_b, 0, size * sizeof(float));
    memset(dens_prev_r, 0, size * sizeof(float));
    memset(dens_prev_g, 0, size * sizeof(float));
    memset(dens_prev_b, 0, size * sizeof(float));

    memset(is_obstacle, 0, size * sizeof(int)); // 默认所有单元格为流体

    return 1;
}

/**
 * @brief 释放流体模拟占用的内存。
 */
static void free_data(void) {
    if (u) free(u);
    if (v) free(v);
    if (u_prev) free(u_prev);
    if (v_prev) free(v_prev);
    if (dens_r) free(dens_r);
    if (dens_g) free(dens_g);
    if (dens_b) free(dens_b);
    if (dens_prev_r) free(dens_prev_r);
    if (dens_prev_g) free(dens_prev_g);
    if (dens_prev_b) free(dens_prev_b);
    if (is_obstacle) free(is_obstacle);
}

/**
 * @brief 为某个量（如速度或密度）添加外部源。
 * @param gridSize 网格大小 N。
 * @param x 目标数组 (当前时刻的值)。
 * @param s 源数组 (外部源的值)。
 * @param dt 时间步长。
 *
 * 改进：不向障碍物单元格添加源。
 */
void add_source(int gridSize, float *x, float *s, float dt) {
    int i, size = (gridSize + 2) * (gridSize + 2);
    for (i = 0; i < size; i++) {
        if (!is_obstacle[i]) { // 仅当不是障碍物时才添加源
            x[i] += dt * s[i];
        }
    }
}

/**
 * @brief 设置边界条件 (仅处理外部边界)。
 * @param gridSize 网格大小 N。
 * @param b 边界类型标识。b=0 表示密度/压力，b=1 表示 u 分量速度，b=2 表示 v 分量速度。
 * @param x 要设置边界的数组。
 */
void set_bnd(int gridSize, int b, float *x) {
    int i;
    for (i = 1; i <= gridSize; i++) {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(gridSize + 1, i)] = b == 1 ? -x[IX(gridSize, i)] : x[IX(gridSize, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, gridSize + 1)] = b == 2 ? -x[IX(i, gridSize)] : x[IX(i, gridSize)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, gridSize + 1)] = 0.5f * (x[IX(1, gridSize + 1)] + x[IX(0, gridSize)]);
    x[IX(gridSize + 1, 0)] = 0.5f * (x[IX(gridSize, 0)] + x[IX(gridSize + 1, 1)]);
    x[IX(gridSize + 1, gridSize + 1)] = 0.5f * (x[IX(gridSize, gridSize + 1)] + x[IX(gridSize + 1, gridSize)]);
}

/**
 * @brief 强制障碍物单元格的数值为0。
 * @param gridSize 网格大小 N。
 * @param x 要处理的数组。
 * @param b 边界类型标识 (0: 密度/压力, 1: u速度, 2: v速度)。
 */
static void enforce_obstacle_boundaries(int gridSize, float *x, int b) {
    int size = (gridSize + 2) * (gridSize + 2);
    for (int i = 0; i < size; i++) {
        if (is_obstacle[i]) {
            x[i] = 0.0f; // 障碍物内部的值设为0
        }
    }
}


/**
 * @brief 泊松方程的通用迭代求解器。
 * 改进：在每次迭代后强制障碍物单元格的值为0。
 */
void lin_solve(int gridSize, int b, float *x, float *x0, float a, float c) {
    for (int k = 0; k < 20; k++) { // 迭代次数
        FOR_EACH_CELL
            if (is_obstacle[IX(i, j)]) continue; // 障碍物单元格不参与计算

            x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                               x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
        END_FOR
        enforce_obstacle_boundaries(gridSize, x, b); // 强制障碍物单元格为0
        set_bnd(gridSize, b, x); // 更新外部边界条件
    }
}

/**
 * @brief 计算流体的扩散项。
 */
void diffuse(int gridSize, int b, float *x, float *x0, float diff, float dt) {
    float a = dt * diff * gridSize * gridSize;
    lin_solve(gridSize, b, x, x0, a, 1 + 4 * a);
}

/**
 * @brief 计算流体的平流项 (Advection)。
 * 改进：在插值时考虑障碍物，并强制障碍物单元格的密度为0。
 */
void advect(int gridSize, int b, float *d, float *d0, float *u, float *v, float dt) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * gridSize;

    FOR_EACH_CELL
        if (is_obstacle[IX(i, j)]) { // 障碍物单元格不计算平流
            d[IX(i, j)] = 0;
            continue;
        }

        x = i - dt0 * u[IX(i, j)];
        y = j - dt0 * v[IX(i, j)];

        if (x < 0.5f) x = 0.5f;
        if (x > gridSize + 0.5f) x = gridSize + 0.5f;
        i0 = (int)x;
        i1 = i0 + 1;

        if (y < 0.5f) y = 0.5f;
        if (y > gridSize + 0.5f) y = gridSize + 0.5f;
        j0 = (int)y;
        j1 = j0 + 1;

        s1 = x - i0;
        s0 = 1 - s1;
        t1 = y - j0;
        t0 = 1 - t1;

        // 双线性插值时，如果插值点落在障碍物内，则将其视为0
        float d00 = is_obstacle[IX(i0, j0)] ? 0.0f : d0[IX(i0, j0)];
        float d10 = is_obstacle[IX(i1, j0)] ? 0.0f : d0[IX(i1, j0)];
        float d01 = is_obstacle[IX(i0, j1)] ? 0.0f : d0[IX(i0, j1)];
        float d11 = is_obstacle[IX(i1, j1)] ? 0.0f : d0[IX(i1, j1)];

        d[IX(i, j)] = s0 * (t0 * d00 + t1 * d01) +
                      s1 * (t0 * d10 + t1 * d11);
    END_FOR

    enforce_obstacle_boundaries(gridSize, d, b); // 强制障碍物单元格为0
    set_bnd(gridSize, b, d); // 更新外部边界条件
}

/**
 * @brief 投射操作 (Project)，用于确保速度场的无散度性（不可压缩性）。
 * 改进：在计算散度时考虑障碍物，并在计算后强制障碍物单元格的速度为0。
 */
void project(int gridSize, float *u, float *v, float *p, float *div) {
    int i, j;

    FOR_EACH_CELL
        if (is_obstacle[IX(i, j)]) { // 障碍物单元格不计算散度
            div[IX(i, j)] = 0;
            p[IX(i, j)] = 0;
            continue;
        }
        
        // 计算散度时，如果相邻单元格是障碍物，其速度视为0
        float u_east = is_obstacle[IX(i + 1, j)] ? 0.0f : u[IX(i + 1, j)];
        float u_west = is_obstacle[IX(i - 1, j)] ? 0.0f : u[IX(i - 1, j)];
        float v_north = is_obstacle[IX(i, j + 1)] ? 0.0f : v[IX(i, j + 1)];
        float v_south = is_obstacle[IX(i, j - 1)] ? 0.0f : v[IX(i, j - 1)];

        div[IX(i, j)] = -0.5f * (u_east - u_west + v_north - v_south) / gridSize;
        p[IX(i, j)] = 0;
    END_FOR

    enforce_obstacle_boundaries(gridSize, div, 0); // 强制障碍物单元格的散度为0
    enforce_obstacle_boundaries(gridSize, p, 0);   // 强制障碍物单元格的压力为0
    set_bnd(gridSize, 0, div);
    set_bnd(gridSize, 0, p);

    // 求解泊松压力方程
    lin_solve(gridSize, 0, p, div, 1, 4);

    FOR_EACH_CELL
        if (is_obstacle[IX(i, j)]) { // 障碍物单元格的速度保持为0
            u[IX(i, j)] = 0;
            v[IX(i, j)] = 0;
            continue;
        }

        // 减去压力梯度时，如果相邻单元格是障碍物，其压力视为0
        float p_east = is_obstacle[IX(i + 1, j)] ? 0.0f : p[IX(i + 1, j)];
        float p_west = is_obstacle[IX(i - 1, j)] ? 0.0f : p[IX(i - 1, j)];
        float p_north = is_obstacle[IX(i, j + 1)] ? 0.0f : p[IX(i, j + 1)];
        float p_south = is_obstacle[IX(i, j - 1)] ? 0.0f : p[IX(i, j - 1)];

        u[IX(i, j)] -= 0.5f * gridSize * (p_east - p_west);
        v[IX(i, j)] -= 0.5f * gridSize * (p_north - p_south);
    END_FOR

    enforce_obstacle_boundaries(gridSize, u, 1); // 强制障碍物单元格的 u 速度为0
    enforce_obstacle_boundaries(gridSize, v, 2); // 强制障碍物单元格的 v 速度为0
    set_bnd(gridSize, 1, u);
    set_bnd(gridSize, 2, v);
}

/**
 * @brief 更新速度场。
 */
void vel_step(int gridSize, float *u, float *v, float *u0, float *v0, float visc, float dt) {
    add_source(gridSize, u, u0, dt);
    add_source(gridSize, v, v0, dt);
    
    SWAP(u0, u);
    diffuse(gridSize, 1, u, u0, visc, dt);

    SWAP(v0, v);
    diffuse(gridSize, 2, v, v0, visc, dt);

    project(gridSize, u, v, u0, v0); // u0, v0 作为临时缓冲区 p 和 div

    SWAP(u0, u);
    SWAP(v0, v);
    advect(gridSize, 1, u, u0, u0, v0, dt);
    advect(gridSize, 2, v, v0, u0, v0, dt);

    project(gridSize, u, v, u0, v0);
}

/**
 * @brief 更新密度值（例如烟雾浓度）。
 */
void dens_step(int gridSize, float *x, float *x0, float *u, float *v, float diff, float dt) {
    add_source(gridSize, x, x0, dt);

    SWAP(x0, x);
    diffuse(gridSize, 0, x, x0, diff, dt);

    SWAP(x0, x);
    advect(gridSize, 0, x, x0, u, v, dt);
}

//--------------------------------------------------------------------------------------
// 障碍物初始化函数
//--------------------------------------------------------------------------------------
static void init_obstacles(void) {
    // 障碍物：水平长条
    int bar_y = N / 2 + 5; // 长条的 y 坐标
    int bar_start_x = N / 4;
    int bar_end_x = 3 * N / 4;
    for (int i = bar_start_x; i <= bar_end_x; i++) {
        is_obstacle[IX(i, bar_y)] = 1;
        is_obstacle[IX(i, bar_y + 1)] = 1; // 增加厚度
        is_obstacle[IX(i, bar_y - 1)] = 1; 
    }

    // 障碍物：圆形 
    int circle_center_x = N / 2;
    int circle_center_y = N / 4; // 圆心位置
    int radius = N / 8;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            float dist_sq = (i - circle_center_x) * (i - circle_center_x) +
                            (j - circle_center_y) * (j - circle_center_y);
            if (dist_sq < radius * radius) {
                is_obstacle[IX(i, j)] = 1;
            }
        }
    }
}


//--------------------------------------------------------------------------------------
// GLUT 相关函数和模拟主循环
//--------------------------------------------------------------------------------------

// 初始化模拟状态
void init_simulation() {
    if (!allocate_data()) {
        exit(1); // 内存分配失败，退出程序
    }
    init_obstacles(); // 初始化障碍物
    printf("流体模拟初始化完成，网格大小: %dx%d\n", N, N);
    printf("--- 键盘控制 ---\n");
    printf("Esc: 退出程序\n");
    printf("R: 设置绘制颜色为红色\n");
    printf("G: 设置绘制颜色为绿色\n");
    printf("B: 设置绘制颜色为蓝色\n");
    printf("W: 设置绘制颜色为白色 (默认)\n");
    printf("Y: 设置绘制颜色为黄色\n");
    printf("C: 设置绘制颜色为青色\n");
    printf("M: 设置绘制颜色为品红色\n");
    printf("--- 鼠标控制 ---\n");
    printf("鼠标左键拖拽: 添加密度源 (当前颜色)\n");
    printf("鼠标右键拖拽: 添加速度源\n");
}

// 模拟步进函数
void simulate_step() {
    int size = (N + 2) * (N + 2);
    // 清空前一时刻的速度和所有颜色通道的密度源
    memset(u_prev, 0, size * sizeof(float));
    memset(v_prev, 0, size * sizeof(float));
    memset(dens_prev_r, 0, size * sizeof(float));
    memset(dens_prev_g, 0, size * sizeof(float));
    memset(dens_prev_b, 0, size * sizeof(float));

    // 用户输入处理（鼠标拖拽）
    // 鼠标左键添加密度，鼠标右键添加速度
    if (mouse_down[0] || mouse_down[2]) {
        // 将窗口坐标转换为网格坐标
        int i = (int)((last_x / (float)win_width) * N) + 1;
        int j = (int)(((win_height - last_y) / (float)win_height) * N) + 1;

        // 确保网格索引在有效范围内，并且不是障碍物
        if (i >= 1 && i <= N && j >= 1 && j <= N && !is_obstacle[IX(i,j)]) {
             if (mouse_down[0]) { // 左键添加密度
                 // 根据当前绘制颜色添加密度源
                 dens_prev_r[IX(i, j)] = current_draw_color[0] * 100.0f; 
                 dens_prev_g[IX(i, j)] = current_draw_color[1] * 100.0f;
                 dens_prev_b[IX(i, j)] = current_draw_color[2] * 100.0f;
             }
             if (mouse_down[2]) { // 右键添加速度
                 // 简单的速度添加，可以根据鼠标移动方向和速度进行改进
                 // 乘以 0.005 是一个经验值，用于控制速度强度
                 u_prev[IX(i, j)] = (last_x - win_width / 2.0) * 0.005; 
                 v_prev[IX(i, j)] = (win_height / 2.0 - last_y) * 0.005; 
             }
        }
    }

    // 更新速度场
    vel_step(N, u, v, u_prev, v_prev, visc, dt);

    // 更新密度场（每个颜色通道独立更新）
    dens_step(N, dens_r, dens_prev_r, u, v, diff, dt);
    dens_step(N, dens_g, dens_prev_g, u, v, diff, dt);
    dens_step(N, dens_b, dens_prev_b, u, v, diff, dt);
}


// GLUT 显示函数
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT); // 清除颜色缓冲区

    glBegin(GL_QUADS); // 绘制四边形
    float h = 1.0f / N; // 每个网格的尺寸 (在 [0, 1] 范围内)

    for (int i = 0; i <= N + 1; i++) {
        for (int j = 0; j <= N + 1; j++) {
            if (is_obstacle[IX(i, j)]) {
                // 如果是障碍物，绘制固定颜色
                glColor3f(0.3f, 0.3f, 0.5f); // 暗蓝灰色
            } else {
                // 否则，使用 RGB 密度值绘制烟雾颜色
                float dr = dens_r[IX(i, j)];
                float dg = dens_g[IX(i, j)];
                float db = dens_b[IX(i, j)];

                // 钳制颜色值在 [0, 1] 范围内
                dr = fmax(0.0f, fmin(1.0f, dr));
                dg = fmax(0.0f, fmin(1.0f, dg));
                db = fmax(0.0f, fmin(1.0f, db));

                glColor3f(dr, dg, db); // 设置颜色 (R, G, B)
            }

            // 绘制当前网格单元的四边形
            // 网格坐标 (i, j) 对应窗口坐标 (i*h, j*h)
            glVertex2f(i * h, j * h);
            glVertex2f((i + 1) * h, j * h);
            glVertex2f((i + 1) * h, (j + 1) * h);
            glVertex2f(i * h, (j + 1) * h);
        }
    }
    glEnd();

    glutSwapBuffers(); // 交换缓冲区，显示绘制结果
}

// GLUT 空闲函数，持续调用以更新模拟和显示
void idle(void) {
    simulate_step(); // 执行一步模拟
    glutPostRedisplay(); // 请求 GLUT 重新绘制窗口
}

// GLUT 窗口大小改变时调用
void reshape(int w, int h) {
    win_width = w;
    win_height = h;
    glViewport(0, 0, w, h); // 设置视口
    glMatrixMode(GL_PROJECTION); // 设置投影矩阵
    glLoadIdentity(); // 加载单位矩阵
    gluOrtho2D(0.0, 1.0, 0.0, 1.0); // 设置二维正交投影，将 [0, 1]x[0, 1] 区域映射到窗口
    glMatrixMode(GL_MODELVIEW); // 设置模型视图矩阵
    glLoadIdentity(); // 加载单位矩阵
}

// GLUT 键盘事件处理
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 27: // 按下 Esc 键退出
            free_data(); // 释放内存
            exit(0);
            break;
        case 'r': case 'R': // 红色
            current_draw_color[0] = 1.0f; current_draw_color[1] = 0.0f; current_draw_color[2] = 0.0f;
            printf("绘制颜色: 红色\n"); break;
        case 'g': case 'G': // 绿色
            current_draw_color[0] = 0.0f; current_draw_color[1] = 1.0f; current_draw_color[2] = 0.0f;
            printf("绘制颜色: 绿色\n"); break;
        case 'b': case 'B': // 蓝色
            current_draw_color[0] = 0.0f; current_draw_color[1] = 0.0f; current_draw_color[2] = 1.0f;
            printf("绘制颜色: 蓝色\n"); break;
        case 'w': case 'W': // 白色 (默认)
            current_draw_color[0] = 1.0f; current_draw_color[1] = 1.0f; current_draw_color[2] = 1.0f;
            printf("绘制颜色: 白色\n"); break;
        case 'y': case 'Y': // 黄色
            current_draw_color[0] = 1.0f; current_draw_color[1] = 1.0f; current_draw_color[2] = 0.0f;
            printf("绘制颜色: 黄色\n"); break;
        case 'c': case 'C': // 青色
            current_draw_color[0] = 0.0f; current_draw_color[1] = 1.0f; current_draw_color[2] = 1.0f;
            printf("绘制颜色: 青色\n"); break;
        case 'm': case 'M': // 品红色
            current_draw_color[0] = 1.0f; current_draw_color[1] = 0.0f; current_draw_color[2] = 1.0f;
            printf("绘制颜色: 品红色\n"); break;
        default: break;
    }
}

// GLUT 鼠标按键事件处理
void mouse(int button, int state, int x, int y) {
    // 记录鼠标按键状态
    if (button < 3) { // 只处理左、中、右键
        mouse_down[button] = (state == GLUT_DOWN);
    }
    last_x = x;
    last_y = y;
}

// GLUT 鼠标移动事件处理 (鼠标按键按下时)
void motion(int x, int y) {
    last_x = x;
    last_y = y;
}


// 主函数
int main(int argc, char **argv) {
    // 初始化 GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); 
    glutInitWindowSize(win_width, win_height); 
    glutCreateWindow("Fliud Simulation"); 

    // 初始化 OpenGL 状态
    glClearColor(0.0, 0.0, 0.0, 1.0); // 设置清除颜色为黑色
    glShadeModel(GL_FLAT); // 设置着色模式为平面着色

    // 初始化模拟数据
    init_simulation();

    // 注册 GLUT 回调函数
    glutDisplayFunc(display); // 注册显示函数
    glutIdleFunc(idle);       // 注册空闲函数
    glutReshapeFunc(reshape); // 注册窗口大小改变函数
    glutKeyboardFunc(keyboard); // 注册键盘事件函数
    glutMouseFunc(mouse);     // 注册鼠标按键事件函数
    glutMotionFunc(motion);   // 注册鼠标移动事件函数

    // 进入 GLUT 主循环
    glutMainLoop();

    return 0;
}