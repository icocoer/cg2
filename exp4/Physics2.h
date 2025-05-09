/**************************************************************************

File: Physics2.h

Prepared by Erkin Tunca for nehe.gamedev.net
(And subsequently modified for keyboard input based on user request)

**************************************************************************/

#ifndef PHYSICS_EXTRA_H // 统一的头文件保护宏 (之前代码中存在多处，已整合)
#define PHYSICS_EXTRA_H

#include "Physics1.h" // Physics1.h is a must for Physics2.h simulations
#include <vector>     // For std::vector
#include <cmath>      // For sqrt, fabs etc.

// class Spring (无修改)
class Spring {
public:
    Mass* mass1;
    Mass* mass2;
    float springConstant;
    float springLength;
    float frictionConstant;

    Spring(Mass* mass1, Mass* mass2,
           float springConstant, float springLength, float frictionConstant) {
        this->springConstant = springConstant;
        this->springLength = springLength;
        this->frictionConstant = frictionConstant;
        this->mass1 = mass1;
        this->mass2 = mass2;
    }

    void solve() {
        if (!mass1 || !mass2) return; // Basic null check

        Vector3D springVector = mass1->pos - mass2->pos;
        float r = springVector.length();
        Vector3D force; // Defaults to (0,0,0)

        if (r != 0) {
            Vector3D e = springVector / r;
            force += e * (r - springLength) * (-springConstant);
            // Dot product for velocity difference along spring axis
            float velDotE = (mass1->vel.x * e.x + mass1->vel.y * e.y + mass1->vel.z * e.z) -
                            (mass2->vel.x * e.x + mass2->vel.y * e.y + mass2->vel.z * e.z);
            force += -e * velDotE * frictionConstant;
        }
        mass1->applyForce(force);
        mass2->applyForce(-force);
    }
};

// class RopeSimulation (无修改)
class RopeSimulation : public Simulation {
public:
    Spring** springs;
    Vector3D gravitation;
    Vector3D ropeConnectionPos;
    Vector3D ropeConnectionVel;

    RopeSimulation(
        int numOfMasses,
        float m,
        float springConstant,
        float springLength,
        float springFrictionConstant,
        Vector3D gravitationVal) : Simulation(numOfMasses, m) {
        this->gravitation = gravitationVal;
        this->ropeConnectionPos = Vector3D(0,0,0); // Initialize
        this->ropeConnectionVel = Vector3D(0,0,0); // Initialize

        if (numOfMasses <= 0) { // Prevent issues with zero or negative masses
             this->springs = nullptr; // Or handle error appropriately
             return;
        }


        for (int a = 0; a < numOfMasses; ++a) {
            if(masses[a]) { // Check if mass exists
                masses[a]->pos.x = a * springLength;
                masses[a]->pos.y = 0;
                masses[a]->pos.z = 0;
            }
        }
        
        if (numOfMasses > 1) {
            springs = new Spring*[numOfMasses - 1];
            for (int a = 0; a < numOfMasses - 1; ++a) {
                 if(masses[a] && masses[a+1]) { // Ensure masses exist before creating spring
                    springs[a] = new Spring(masses[a], masses[a + 1],
                                        springConstant, springLength, springFrictionConstant);
                 } else {
                    springs[a] = nullptr; // Handle missing mass case
                 }
            }
        } else {
            springs = nullptr; // No springs if less than 2 masses
        }
    }

    void release() override {
        if (springs) {
             for (int a = 0; a < numOfMasses - 1; ++a) {
                if (springs[a]) {
                    delete springs[a];
                    springs[a] = NULL;
                }
            }
            delete[] springs; // Use delete[] for arrays
            springs = NULL;
        }
        Simulation::release();
    }

    void solve() override {
        if (springs) {
            for (int a = 0; a < numOfMasses - 1; ++a) {
                if (springs[a]) {
                    springs[a]->solve();
                }
            }
        }

        for (int a = 0; a < numOfMasses; ++a) {
            if(masses[a]) {
                masses[a]->applyForce(gravitation * masses[a]->m);
            }
        }
    }

    void simulate(float dt) override {
        Simulation::simulate(dt);
        ropeConnectionPos += ropeConnectionVel * dt;
        if (numOfMasses > 0 && masses[0]) { // Ensure first mass exists
            masses[0]->pos = ropeConnectionPos;
            masses[0]->vel = ropeConnectionVel;
        }
    }

    void setRopeConnectionVel(Vector3D newRopeConnectionVel) {
        this->ropeConnectionVel = newRopeConnectionVel;
    }
};


// --- ElasticSolidSimulation Class (修改部分) ---
class ElasticSolidSimulation : public Simulation {
public:
    int dimX, dimY, dimZ;
    std::vector<Spring*> springs; // 使用 std::vector 更安全方便
    Vector3D gravitation;
    Vector3D keyboardAppliedForce; // 用于存储键盘施加的力

    float groundLevel;
    float restitutionCoefficient;
    float staticFrictionCoefficient;
    float dynamicFrictionCoefficient;
    float timeStep; // 用于物理计算的时间步长

    ElasticSolidSimulation(
        int dX, int dY, int dZ,
        float massVal,
        float springConstStructural, float springLenStructural,
        float springConstShear, float springLenShearFactor,
        float springConstBend, float springLenBendFactor,
        float frictionConstSpring,
        Vector3D grav,
        float gLevel, float restitution,
        float staticFriction, float dynamicFriction,
        float dt // 从 main.cpp 传入 timeStep
    ) : Simulation(dX * dY * dZ, massVal),
        dimX(dX), dimY(dY), dimZ(dZ), gravitation(grav),
        groundLevel(gLevel), restitutionCoefficient(restitution),
        staticFrictionCoefficient(staticFriction), dynamicFrictionCoefficient(dynamicFriction),
        timeStep(dt) // 初始化 timeStep
    {
        if (numOfMasses == 0) return;

        keyboardAppliedForce = Vector3D(0, 0, 0); // 初始化键盘力

        // 1. Initialize Mass Positions
        float initialHeightOffset = (dimY * springLenStructural * 0.5f) + springLenStructural * 2.0f;
        for (int k = 0; k < dimZ; ++k) {
            for (int j = 0; j < dimY; ++j) {
                for (int i = 0; i < dimX; ++i) {
                    Mass* currentMass = getMassByIndex(i, j, k);
                    if (currentMass) {
                        currentMass->pos = Vector3D(
                            i * springLenStructural - (dimX - 1) * springLenStructural * 0.5f,
                            j * springLenStructural + initialHeightOffset,
                            k * springLenStructural - (dimZ - 1) * springLenStructural * 0.5f
                        );
                        currentMass->vel = Vector3D(0, 0, 0);
                    }
                }
            }
        }

        // 2. Create Springs
        for (int k = 0; k < dimZ; ++k) {
            for (int j = 0; j < dimY; ++j) {
                for (int i = 0; i < dimX; ++i) {
                    Mass* m1 = getMassByIndex(i, j, k);
                    if (!m1) continue;

                    // Structural
                    if (i + 1 < dimX) springs.push_back(new Spring(m1, getMassByIndex(i + 1, j, k), springConstStructural, springLenStructural, frictionConstSpring));
                    if (j + 1 < dimY) springs.push_back(new Spring(m1, getMassByIndex(i, j + 1, k), springConstStructural, springLenStructural, frictionConstSpring));
                    if (k + 1 < dimZ) springs.push_back(new Spring(m1, getMassByIndex(i, j, k + 1), springConstStructural, springLenStructural, frictionConstSpring));

                    // Shear
                    float shearRestLength = springLenStructural * springLenShearFactor;
                    if (springConstShear > 0 && shearRestLength > 0) {
                        if (i + 1 < dimX && j + 1 < dimY) {
                            springs.push_back(new Spring(m1, getMassByIndex(i + 1, j + 1, k), springConstShear, shearRestLength, frictionConstSpring));
                            Mass* m_i1_j0 = getMassByIndex(i + 1, j, k); Mass* m_i0_j1 = getMassByIndex(i, j + 1, k);
                            if (m_i1_j0 && m_i0_j1) springs.push_back(new Spring(m_i1_j0, m_i0_j1, springConstShear, shearRestLength, frictionConstSpring));
                        }
                        if (i + 1 < dimX && k + 1 < dimZ) {
                            springs.push_back(new Spring(m1, getMassByIndex(i + 1, j, k + 1), springConstShear, shearRestLength, frictionConstSpring));
                            Mass* m_i1_j0 = getMassByIndex(i + 1, j, k); Mass* m_i0_k1 = getMassByIndex(i, j, k + 1);
                            if (m_i1_j0 && m_i0_k1) springs.push_back(new Spring(m_i1_j0, m_i0_k1, springConstShear, shearRestLength, frictionConstSpring));
                        }
                        if (j + 1 < dimY && k + 1 < dimZ) {
                            springs.push_back(new Spring(m1, getMassByIndex(i, j + 1, k + 1), springConstShear, shearRestLength, frictionConstSpring));
                            Mass* m_i0_j1 = getMassByIndex(i, j + 1, k); Mass* m_i0_k1 = getMassByIndex(i, j, k + 1);
                            if (m_i0_j1 && m_i0_k1) springs.push_back(new Spring(m_i0_j1, m_i0_k1, springConstShear, shearRestLength, frictionConstSpring));
                        }
                    }
                    // Bending
                    float bendRestLength = springLenStructural * springLenBendFactor;
                     if (springConstBend > 0 && bendRestLength > 0) {
                        if (i + 2 < dimX) springs.push_back(new Spring(m1, getMassByIndex(i + 2, j, k), springConstBend, bendRestLength, frictionConstSpring));
                        if (j + 2 < dimY) springs.push_back(new Spring(m1, getMassByIndex(i, j + 2, k), springConstBend, bendRestLength, frictionConstSpring));
                        if (k + 2 < dimZ) springs.push_back(new Spring(m1, getMassByIndex(i, j, k + 2), springConstBend, bendRestLength, frictionConstSpring));
                    }
                }
            }
        }
    }

    // setExternalForceForKey 方法定义
    void setExternalForceForKey(const Vector3D& force) {
        keyboardAppliedForce = force;
    }

    Mass* getMassByIndex(int x, int y, int z) {
        if (x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ) {
            return nullptr;
        }
        int index = x + y * dimX + z * dimX * dimY;
        return (index >= 0 && index < numOfMasses) ? masses[index] : nullptr;
    }

    void solve() override {
        // 1. 计算弹簧力
        for (Spring* s : springs) {
            if (s) s->solve();
        }

        // 2. 对每个质量点施加其他力
        for (int i = 0; i < numOfMasses; ++i) {
            Mass* currentMass = masses[i];
            if (!currentMass) continue;

            // a. 施加重力
            currentMass->applyForce(gravitation * currentMass->m);

            // b. 施加键盘输入的力
            if (keyboardAppliedForce.x != 0.0f || keyboardAppliedForce.y != 0.0f || keyboardAppliedForce.z != 0.0f) {
                currentMass->applyForce(keyboardAppliedForce);
            }

            // c. 处理地面碰撞和摩擦
            if (currentMass->pos.y < groundLevel + 0.01f) { // 使用一个小的epsilon比较
                currentMass->pos.y = groundLevel; // 修正位置，防止穿透

                float incomingNormalForceY = 0;
                if (currentMass->force.y < 0) { // 如果有力正把它推向地面
                    incomingNormalForceY = -currentMass->force.y; // Y方向的法向力分量
                }


                // 碰撞响应 (速度反弹)
                if (currentMass->vel.y < 0) {
                    currentMass->vel.y *= -restitutionCoefficient;
                }

                // 计算摩擦力
                Vector3D tangentialVelocity = Vector3D(currentMass->vel.x, 0, currentMass->vel.z);
                float tangentialSpeed = tangentialVelocity.length();

                if (tangentialSpeed > 0.0001f) { // 如果有切向速度
                    Vector3D frictionDirection = -tangentialVelocity.unit(); // 摩擦力方向与切向速度相反

                    // 简化的法向力 N 的计算 (重力 + 其他压向地面的力)
                    // 更精确的模型会考虑所有垂直于地面的力的总和
                    float N = currentMass->m * std::abs(gravitation.y) + incomingNormalForceY;
                    if (N <= 0) { // 确保法向力为正，至少有重力支撑
                        N = currentMass->m * std::abs(gravitation.y);
                         if (N <= 0) N = 1.0f; // 避免除零，如果重力也为零或负（不太可能）
                    }


                    float dynamicFrictionMagnitude = dynamicFrictionCoefficient * N;
                    Vector3D frictionForce = frictionDirection * dynamicFrictionMagnitude;

                    // 施加摩擦力前，确保摩擦力不会在一个时间步内反转速度
                    Vector3D velChangeDueToFriction = (frictionForce / currentMass->m) * this->timeStep; // 使用成员变量 timeStep
                    if (velChangeDueToFriction.length() >= tangentialSpeed) {
                        // 施加恰好能使其切向速度变为0的力
                        Vector3D exactFriction = -tangentialVelocity * (currentMass->m / this->timeStep);
                        currentMass->applyForce(exactFriction);
                        currentMass->vel.x = 0; // 速度直接置零可能更稳定
                        currentMass->vel.z = 0;

                    } else {
                        currentMass->applyForce(frictionForce);
                    }
                }
            }
        }
    }

    void release() override {
        for (Spring* s : springs) {
            if (s) delete s;
        }
        springs.clear();
        Simulation::release(); // 调用基类的 release
    }
};

#endif // PHYSICS_EXTRA_H