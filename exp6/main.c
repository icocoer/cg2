#include <stdio.h>
#include <stdlib.h>
#include <string.h> // 用于 memset, memcpy 等内存操作函数
#include <GL/glut.h> // 包含 GLUT 库

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
// dens: 当前时刻的密度值 (例如烟雾浓度)
// dens_prev: 前一时刻的密度值，也用作临时缓冲区
static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;

// 模拟参数
static float dt = 0.1f; // 时间步长 (delta t)
// 粘度系数 (viscosity)，衡量流体抵抗形变的能力。
// 在扩散项计算中使用。
static float visc = 0.0000001f; // 设为一个很小的值以减少扩散影响
// 密度扩散系数 (diffusion)，衡量密度扩散速率。
// 在密度扩散计算中使用。
static float diff = 0.0000001f; // 设为一个很小的值以减少扩散影响

// 窗口尺寸
static int win_width = 512;
static int win_height = 512;

//--------------------------------------------------------------------------------------
// 函数声明
//--------------------------------------------------------------------------------------
// 将参数 N 重命名为 gridSize，避免与宏 N 冲突
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
    dens = (float *)malloc(size * sizeof(float));
    dens_prev = (float *)malloc(size * sizeof(float));

    if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
        fprintf(stderr, "无法分配数据内存\n");
        return 0;
    }

    // 初始化所有数据为0
    memset(u, 0, size * sizeof(float));
    memset(v, 0, size * sizeof(float));
    memset(u_prev, 0, size * sizeof(float));
    memset(v_prev, 0, size * sizeof(float));
    memset(dens, 0, size * sizeof(float));
    memset(dens_prev, 0, size * sizeof(float));

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
    if (dens) free(dens);
    if (dens_prev) free(dens_prev);
}

/**
 * @brief 为某个量（如速度或密度）添加外部源。
 * @param gridSize 网格大小 N。
 * @param x 目标数组 (当前时刻的值)。
 * @param s 源数组 (外部源的值)。
 * @param dt 时间步长。
 */
void add_source(int gridSize, float *x, float *s, float dt) {
    int i, size = (gridSize + 2) * (gridSize + 2);
    for (i = 0; i < size; i++) {
        // 数学运算：x[i] = x[i] + dt * s[i]
        // 解释：在当前时间步长 dt 内，将源 s[i] 乘以 dt 累加到 x[i] 上。
        // 这对应于 Navier-Stokes 方程中的外力项 F 或密度方程中的源项 S，表示外部因素对流体状态的改变。
        x[i] += dt * s[i];
    }
}

/**
 * @brief 设置边界条件。
 * @param gridSize 网格大小 N。
 * @param b 边界类型标识。b=0 表示密度/压力，b=1 表示 u 分量速度，b=2 表示 v 分量速度。
 * @param x 要设置边界的数组。
 *
 * 根据 Stam 的实现，当 b=0 时，对密度和压力场应用 Neumann 边界条件（法向导数为 0），
 * 即边界值等于相邻内部网格的值。
 * 当 b=1 或 b=2 时，对速度场应用无滑移 (no-slip) 边界条件，即边界速度为 0。
 * 具体实现中，通过设置边界处的速度为相邻内部网格速度的负值来实现无滑移。
 * 角点的处理是取相邻两个边界网格的平均值。
 */
void set_bnd(int gridSize, int b, float *x) {
    int i;
    for (i = 1; i <= gridSize; i++) {
        // 左右边界 (x 方向速度 u，法向导数为 0 的密度/压力)
        // x[IX(0, i)]: 左边界，x[IX(1, i)]: 左边界内侧第一个网格
        // 当 b=1 (u 速度) 时，x[IX(0,i)] = -x[IX(1,i)] 实现 u=0 (无滑移条件)。
        // 当 b!=1 (v 速度, 密度, 压力) 时，x[IX(0,i)] = x[IX(1,i)] 实现法向导数为 0 (Neumann 条件)。
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(gridSize + 1, i)] = b == 1 ? -x[IX(gridSize, i)] : x[IX(gridSize, i)];

        // 上下边界 (y 方向速度 v，法向导数为 0 的密度/压力)
        // x[IX(i, 0)]: 下边界，x[IX(i, 1)]: 下边界内侧第一个网格
        // 当 b=2 (v 速度) 时，x[IX(i,0)] = -x[IX(i,1)] 实现 v=0 (无滑移条件)。
        // 当 b!=2 (u 速度, 密度, 压力) 时，x[IX(i,0)] = x[IX(i,1)] 实现法向导数为 0 (Neumann 条件)。
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, gridSize + 1)] = b == 2 ? -x[IX(i, gridSize)] : x[IX(i, gridSize)];
    }

    // 角点处理
    // 角点值取其相邻两个边界点的平均值，这有助于在迭代过程中稳定边界条件。
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, gridSize + 1)] = 0.5f * (x[IX(1, gridSize + 1)] + x[IX(0, gridSize)]);
    x[IX(gridSize + 1, 0)] = 0.5f * (x[IX(gridSize, 0)] + x[IX(gridSize + 1, 1)]);
    x[IX(gridSize + 1, gridSize + 1)] = 0.5f * (x[IX(gridSize, gridSize + 1)] + x[IX(gridSize + 1, gridSize)]);
}

/**
 * @brief 泊松方程的通用迭代求解器。
 * @param gridSize 网格大小 N。
 * @param b 边界类型标识。
 * @param x 待求解的数组 (当前迭代值)。
 * @param x0 右侧项/源项数组 (初始值或上一时间步的值)。
 * @param a 迭代公式中的参数 a。
 * @param c 迭代公式中的参数 c。
 *
 * 该函数用于求解形如 Ax=b 的线性方程组，特别是扩散项和压力泊松方程。
 * 迭代公式来自泊松方程和黏度扩散方程的离散化形式。
 * 通用形式：$x_{i,j}^{(k+1)} = \frac{(x_{i+1,j}^{(k)} + x_{i-1,j}^{(k)} + x_{i,j+1}^{(k)} + x_{i,j-1}^{(k)}) + \alpha b_{i,j}}{\beta}$。
 * 这里 `lin_solve` 内部的迭代公式为：
 * $x_{i,j}^{(k+1)} = \frac{x0_{i,j} + a \times (x_{i-1,j}^{(k)} + x_{i+1,j}^{(k)} + x_{i,j-1}^{(k)} + x_{i,j+1}^{(k)})}{c}$。
 * 它通过雅可比迭代或高斯-赛德尔迭代（此处由于循环顺序，更像高斯-赛德尔）进行求解。
 * 迭代 20 次是为了收敛到近似解。
 */
void lin_solve(int gridSize, int b, float *x, float *x0, float a, float c) {
    for (int k = 0; k < 20; k++) { // 迭代次数，通常 20 次足以获得稳定结果
        FOR_EACH_CELL
            // 数学运算：x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
            // 解释：求解当前网格 (i,j) 的值，基于其旧值 x0[IX(i,j)] 和相邻网格的当前/旧值。
            // 这是泊松方程 $\nabla^2 p = \text{div}$ 或黏度扩散方程 $(I - \nu \Delta t \nabla^2) u = u_{prev}$ 离散化后通过迭代方法求解的公式。
            // 对于扩散项：x0 对应于上一时刻的速度 u_prev，a = $\nu \Delta t / (\Delta x)^2$，c = 1 + 4a。
            // 对于压力泊松方程：x0 对应于散度 div，a = $(\Delta x)^2$，c = 4。
            // Stam 的简化实现中，通过调整传入的 a 和 c 来使用同一个函数求解这两种方程。
            x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                               x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
        END_FOR
        set_bnd(gridSize, b, x); // 每次迭代后更新边界条件，确保边界值符合要求
    }
}

/**
 * @brief 计算流体的扩散项。
 * @param gridSize 网格大小 N。
 * @param b 边界类型标识。
 * @param x 目标数组 (当前时刻的值)。
 * @param x0 源数组 (上一时刻的值)。
 * @param diff 扩散系数 (可以是粘度 visc 或密度扩散系数 diff)。
 * @param dt 时间步长。
 *
 * 扩散项的偏微分方程为：$\frac{\partial q}{\partial t} = \kappa \nabla^{2} q$，其中 q 是待扩散的量，$\kappa$ 是扩散系数。
 * Stam 采用隐式积分形式：$(I - \kappa \Delta t \nabla^{2}) q(x, t+\Delta t) = q(x, t)$。
 * 离散化后，这是一个泊松方程，通过 lin_solve 函数求解。
 */
void diffuse(int gridSize, int b, float *x, float *x0, float diff, float dt) {
    // 数学运算：a = dt * diff * gridSize * gridSize;
    // 解释：计算迭代求解器中的参数 a。
    // 在扩散方程离散化中，$(\Delta x)^2 = (1/gridSize)^2 = 1/gridSize^2$。
    // 迭代公式中的 a 对应 $\kappa \Delta t / (\Delta x)^2 = \text{diff} \times \text{dt} \times \text{gridSize}^2$。
    float a = dt * diff * gridSize * gridSize;
    // 调用通用泊松求解器 lin_solve。
    // c 的值为 1 + 4*a，对应于扩散方程离散化后系数矩阵对角线元素。
    // 迭代公式为：$x_{i,j}^{(k+1)} = \frac{x0_{i,j} + a \times (x_{i-1,j}^{(k)} + x_{i+1,j}^{(k)} + x_{i,j-1}^{(k)} + x_{i,j+1}^{(k)})}{1 + 4a}$。
    lin_solve(gridSize, b, x, x0, a, 1 + 4 * a);
}

/**
 * @brief 计算流体的平流项 (Advection)。
 * @param gridSize 网格大小 N。
 * @param b 边界类型标识。
 * @param d 目标数组 (当前时刻的值)。
 * @param d0 源数组 (上一时刻的值)。
 * @param u 速度场 u 分量。
 * @param v 速度场 v 分量。
 * @param dt 时间步长。
 *
 * 平流项描述了流体中某个量（如速度或密度）如何沿着速度场移动。
 * Stam 采用隐式积分法 (Backwards Tracing)，从当前点 P 逆向追踪粒子轨道，找到它在上一时刻 P' 的位置，
 * 然后将 P' 处的量（通过双线性插值获得）复制到 P 点。
 * 数学公式：$q(x, t+\Delta t) = q(x - u(x, t) \Delta t, t)$，其中 q 是待平流的量，u 是速度场。
 */
void advect(int gridSize, int b, float *d, float *d0, float *u, float *v, float dt) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    // 数学运算：dt0 = dt * gridSize;
    // 解释：将时间步长 dt 乘以 gridSize，将物理距离转换为网格距离。
    // 速度场 u, v 通常表示每单位物理距离的速度，乘以 dt 得到物理位移，再乘以 gridSize 得到网格位移。
    dt0 = dt * gridSize;

    FOR_EACH_CELL
        // 数学运算：x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
        // 解释：逆向追踪粒子位置。从当前网格中心 (i,j) 出发，沿着速度场 (u,v) 的负方向移动 dt0 的距离，找到粒子在上一时刻的源点 (x, y)。
        // 这里 (i,j) 是网格中心在整数网格坐标系中的位置。
        x = i - dt0 * u[IX(i, j)];
        y = j - dt0 * v[IX(i, j)];

        // 边界钳制：确保追踪到的源点位置 (x, y) 在有效网格范围内 [0.5, gridSize+0.5]。
        // 0.5f 和 gridSize+0.5f 是为了避免整数截断和确保插值点总是在网格中心之间。
        if (x < 0.5f) x = 0.5f;
        if (x > gridSize + 0.5f) x = gridSize + 0.5f;
        i0 = (int)x; // 源点 x 坐标的整数部分，作为左下角网格的 x 索引。
        i1 = i0 + 1; // 右上角网格的 x 索引。

        if (y < 0.5f) y = 0.5f;
        if (y > gridSize + 0.5f) y = gridSize + 0.5f;
        j0 = (int)y; // 源点 y 坐标的整数部分，作为左下角网格的 y 索引。
        j1 = j0 + 1; // 右上角网格的 y 索引。

        // 双线性插值权重计算。
        // s1 是 x 在网格单元内的相对位置 (小数部分)，s0 = 1 - s1 是其补数。
        // t1 是 y 在网格单元内的相对位置 (小数部分)，t0 = 1 - t1 是其补数。
        s1 = x - i0;
        s0 = 1 - s1;
        t1 = y - j0;
        t0 = 1 - t1;

        // 数学运算：双线性插值
        // d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        // 解释：根据源点 (x, y) 周围四个相邻网格点 (i0,j0), (i1,j0), (i0,j1), (i1,j1) 在上一时刻的值 d0 进行双线性插值。
        // 首先在 x 方向进行两次线性插值：
        // 在 j0 行，(x, j0) 处的值近似为：$d_{x,j0} = s0 \times d0_{i0,j0} + s1 \times d0_{i1,j0}$
        // 在 j1 行，(x, j1) 处的值近似为：$d_{x,j1} = s0 \times d0_{i0,j1} + s1 \times d0_{i1,j1}$
        // 然后在 y 方向进行一次线性插值：
        // 在 x 列，(x, y) 处的值近似为：$d_{x,y} = t0 \times d_{x,j0} + t1 \times d_{x,j1}$
        // 将上面两步合并，得到最终的双线性插值公式。
        d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                      s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    END_FOR

    set_bnd(gridSize, b, d); // 更新边界条件
}

/**
 * @brief 投射操作 (Project)，用于确保速度场的无散度性（不可压缩性）。
 * @param gridSize 网格大小 N。
 * @param u 速度场 u 分量 (待修正)。
 * @param v 速度场 v 分量 (待修正)。
 * @param p 压力场 (临时数组，用于存储泊松方程的解)。
 * @param div 散度场 (临时数组，用于存储速度场的散度)。
 *
 * 根据 Helmholtz-Hodage 分解定律，任何矢量场 w 都可以分解为一个无散度矢量场 u
 * 和一个标量场 p 的梯度之和：w = u + $\nabla$p。
 * 为了得到无散度的速度场 u，我们可以计算 u = w - $\nabla$p。
 * 其中，压力 p 可以通过求解泊松方程 $\nabla^{2}p = \nabla \cdot w$ 得到。
 * 这里的 w 是经过平流、扩散和外力作用后的临时速度场。
 */
void project(int gridSize, float *u, float *v, float *p, float *div) {
    int i, j;

    FOR_EACH_CELL
        // 数学运算：div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / gridSize;
        // 解释：计算速度场 w 的散度 (divergence) $\nabla \cdot w = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}$。
        // 采用中心差分法近似导数。$\frac{\partial u}{\partial x}$ 近似为 $(u_{i+1,j} - u_{i-1,j}) / (2\Delta x)$。
        // 这里的 $\Delta x = 1/gridSize$，所以 $1/(2\Delta x) = gridSize/2$。
        // 公式中的 -0.5f/gridSize 来自于泊松方程 $\nabla^{2}p = \nabla \cdot w$ 中的项。
        // 这里的 div 数组存储的是 $(\Delta x \Delta y) \nabla \cdot w$ 的近似值。
        div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] +
                                 v[IX(i, j + 1)] - v[IX(i, j - 1)]) / gridSize;
        p[IX(i, j)] = 0; // 初始化压力场 p 为 0。
    END_FOR

    // 设置散度场和压力场的边界条件。
    set_bnd(gridSize, 0, div);
    set_bnd(gridSize, 0, p);

    // 求解泊松压力方程 $\nabla^{2}p = \nabla \cdot w$。
    // 调用 lin_solve 求解压力 p。
    // 离散化形式为：$\frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{(\Delta x)^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{(\Delta y)^2} = \nabla \cdot w_{i,j}$
    // 若 $\Delta x = \Delta y$，则 $\frac{p_{i+1,j} + p_{i-1,j} + p_{i,j+1} + p_{i,j-1} - 4p_{i,j}}{(\Delta x)^2} = \nabla \cdot w_{i,j}$
    // 整理得：$4p_{i,j} - (p_{i+1,j} + p_{i-1,j} + p_{i,j+1} + p_{i,j-1}) = -(\Delta x)^2 \nabla \cdot w_{i,j}$
    // 迭代公式：$p_{i,j}^{(k+1)} = \frac{(p_{i+1,j}^{(k)} + p_{i-1,j}^{(k)} + p_{i,j+1}^{(k)} + p_{i,j-1}^{(k)}) - (\Delta x)^{2}\nabla\cdot w_{i,j}}{4}$。
    // Stam 的实现中，lin_solve(gridSize, 0, p, div, 1, 4) 对应于迭代公式中的 x=p, x0=div, a=1, c=4。
    // 这里的 div 存储的是 $(\Delta x \Delta y) \nabla \cdot w$ 的近似值，所以迭代公式略有不同，但数值上是等价的。
    lin_solve(gridSize, 0, p, div, 1, 4);

    FOR_EACH_CELL
        // 数学运算：u[IX(i, j)] -= 0.5f * gridSize * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
        // 数学运算：v[IX(i, j)] -= 0.5f * gridSize * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
        // 解释：根据 Helmholtz-Hodage 分解，从原始速度场 w 中减去压力场的梯度 $\nabla p$ 得到无散度速度场 u。
        // u_new = w - $\nabla p$。 $\nabla p = (\frac{\partial p}{\partial x}, \frac{\partial p}{\partial y})$。
        // $\frac{\partial p}{\partial x}$ 近似为 $(p_{i+1,j} - p_{i-1,j}) / (2\Delta x)$。
        // 这里的 $0.5f * gridSize$ 对应 $1/(2\Delta x)$。
        // w 在这里就是经过前面步骤（外力、扩散、平流）后的 u 和 v 数组。
        u[IX(i, j)] -= 0.5f * gridSize * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
        v[IX(i, j)] -= 0.5f * gridSize * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    END_FOR

    // 修正后的速度场 u, v 再次设置边界条件，确保无散度速度场也满足边界条件。
    set_bnd(gridSize, 1, u);
    set_bnd(gridSize, 2, v);
}

/**
 * @brief 更新速度场。
 * @param gridSize 网格大小 N。
 * @param u 当前时刻 u 分量。
 * @param v 当前时刻 v 分量。
 * @param u0 前一时刻 u 分量 (用作输入源和临时缓冲区)。
 * @param v0 前一时刻 v 分量 (用作输入源和临时缓冲区)。
 * @param visc 粘度系数。
 * @param dt 时间步长。
 *
 * 速度场更新的步骤 (根据 Navier-Stokes 方程分解)：
 * 1. 添加外力 (Add Forces)。
 * 2. 计算粘性扩散 (Viscous Diffuse)。
 * 3. 投射操作 (Project) 消除散度。
 * 4. 计算自身平流 (Self-advection)。
 * 5. 再次投射操作 (Project) 消除散度，确保守恒性。
 */
void vel_step(int gridSize, float *u, float *v, float *u0, float *v0, float visc, float dt) {
    // 1. 添加外力：将前一时刻的速度场 u0, v0 作为外部源添加到当前速度场 u, v 中。
    // 注意：这里的 u0, v0 实际上是用户输入或模拟中的外部力/速度增量。
    add_source(gridSize, u, u0, dt);
    add_source(gridSize, v, v0, dt);

    // 2. 粘性扩散：计算速度场的扩散。
    // 先交换 u 和 u0 的指针，使得 u0 存储当前速度，u 作为计算扩散结果的临时缓冲区。
    SWAP(u0, u);
    diffuse(gridSize, 1, u, u0, visc, dt); // 对 u 分量进行扩散计算，边界类型 b=1

    SWAP(v0, v);
    diffuse(gridSize, 2, v, v0, visc, dt); // 对 v 分量进行扩散计算，边界类型 b=2

    // 3. 投射操作：消除速度场中的散度，确保不可压缩性。
    // u0, v0 在这里作为 project 函数内部临时存储压力 p 和散度 div 的缓冲区。
    project(gridSize, u, v, u0, v0);

    // 4. 自身平流：计算速度场随自身流动的平流效应。
    // 再次交换 u 和 u0 的指针，使得 u0 存储扩散后的速度，u 作为计算平流结果的临时缓冲区。
    // u0, v0 (交换后) 现在是平流运算的输入速度场 (advecting field)。
    SWAP(u0, u);
    SWAP(v0, v);
    advect(gridSize, 1, u, u0, u0, v0, dt); // 对 u 分量进行平流计算，边界类型 b=1
    advect(gridSize, 2, v, v0, u0, v0, dt); // 对 v 分量进行平流计算，边界类型 b=2

    // 5. 再次投射操作：Stam 建议在平流之后再次进行投射，以确保数值稳定性并维持散度为零。
    project(gridSize, u, v, u0, v0);
}

/**
 * @brief 更新密度值（例如烟雾浓度）。
 * @param gridSize 网格大小 N。
 * @param x 当前时刻密度值。
 * @param x0 前一时刻密度值 (用作输入源和临时缓冲区)。
 * @param u 速度场 u 分量。
 * @param v 速度场 v 分量。
 * @param diff 密度扩散系数。
 * @param dt 时间步长。
 *
 * 密度更新的步骤 (根据密度更新方程)：
 * 1. 添加密度源 (Add Sources)。
 * 2. 计算密度扩散 (Diffuse)。
 * 3. 计算密度平流 (Move)。
 */
void dens_step(int gridSize, float *x, float *x0, float *u, float *v, float diff, float dt) {
    // 1. 添加密度源：将前一时刻的密度源 x0 累加到当前密度 x 中。
    add_source(gridSize, x, x0, dt);

    // 2. 密度扩散：计算密度场的扩散。
    SWAP(x0, x);
    diffuse(gridSize, 0, x, x0, diff, dt); // 对密度进行扩散计算，边界类型 b=0

    // 3. 密度平流：计算密度场随速度场 u, v 的平流效应。
    SWAP(x0, x);
    advect(gridSize, 0, x, x0, u, v, dt); // 对密度进行平流计算，边界类型 b=0
}


//--------------------------------------------------------------------------------------
// GLUT 相关函数和模拟主循环
//--------------------------------------------------------------------------------------

// 用于存储鼠标拖拽信息
static int mouse_down[3]; // 鼠标按键状态
static int last_x, last_y; // 上一帧鼠标位置

// 初始化模拟状态
void init_simulation() {
    if (!allocate_data()) {
        exit(1); // 内存分配失败，退出程序
    }
    printf("流体模拟初始化完成，网格大小: %dx%d\n", N, N);
}

// 模拟步进函数
void simulate_step() {
    // 在这里处理用户输入（鼠标拖拽）作为速度和密度源
    int size = (N + 2) * (N + 2);
    memset(u_prev, 0, size * sizeof(float));
    memset(v_prev, 0, size * sizeof(float));
    memset(dens_prev, 0, size * sizeof(float));

    // 示例：如果在窗口中央添加一个持续的密度源
    // dens_prev[IX(N/2, N/2)] = 100.0f;

    // 用户输入处理（鼠标拖拽）
    // 鼠标左键添加密度，鼠标右键添加速度
    if (mouse_down[0] || mouse_down[2]) {
        // 将窗口坐标转换为网格坐标
        int i = (int)((last_x / (float)win_width) * N) + 1;
        int j = (int)(((win_height - last_y) / (float)win_height) * N) + 1;

        // 确保网格索引在有效范围内
        if (i >= 1 && i <= N && j >= 1 && j <= N) {
             if (mouse_down[0]) { // 左键添加密度
                 dens_prev[IX(i, j)] = 100.0f; // 添加密度源
             }
             if (mouse_down[2]) { // 右键添加速度
                 // 简单的速度添加，可以根据鼠标移动方向和速度进行改进
                 u_prev[IX(i, j)] = (last_x - win_width / 2.0) * 0.005; // 示例速度源
                 v_prev[IX(i, j)] = (win_height / 2.0 - last_y) * 0.005; // 示例速度源
             }
        }
    }


    // 更新速度场
    vel_step(N, u, v, u_prev, v_prev, visc, dt);

    // 更新密度场
    dens_step(N, dens, dens_prev, u, v, diff, dt);
}


// GLUT 显示函数
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT); // 清除颜色缓冲区

    // 使用密度值绘制网格
    glBegin(GL_QUADS); // 绘制四边形
    float h = 1.0f / N; // 每个网格的尺寸 (在 [0, 1] 范围内)

    for (int i = 0; i <= N + 1; i++) {
        for (int j = 0; j <= N + 1; j++) {
            // 获取当前网格的密度值，并钳制在 [0, 1] 范围内用于颜色映射
            float d = dens[IX(i, j)];
            if (d < 0.0f) d = 0.0f;
            if (d > 1.0f) d = 1.0f; // 可以根据需要调整密度到颜色的映射范围

            // 将密度映射到灰度（或其他颜色）
            // 这里使用灰度，密度越高颜色越亮
            glColor3f(d, d, d); // 设置颜色 (R, G, B)

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
    if (key == 27) { // 按下 Esc 键退出
        free_data(); // 释放内存
        exit(0);
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
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // 双缓冲区，RGB 颜色模式
    glutInitWindowSize(win_width, win_height); // 设置窗口尺寸
    glutCreateWindow("流体模拟"); // 创建窗口并设置标题

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
