/*
Computer animation assignment 3: Inverse Kinematics

First step: Search TODO comments to find the methods that you need to implement.
    - src/simulation/kinematics.cpp
*/
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>

#include "Eigen/Core"
#include "GLFW/glfw3.h"
#define GLAD_GL_IMPLEMENTATION
#include "glad/gl.h"
#undef GLAD_GL_IMPLEMENTATION
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "acclaim.h"
#include "graphics.h"
#include "simulation.h"
#include "util.h"

// Global variables are evil!
namespace {
// Shadow texture size, default is 2048 * 2048
int shadowTextureSize = 2048;
// Scene size
int g_ScreenWidth = 1024, g_ScreenHeight = 768;
// Camera controlled by panel
graphics::DefaultCamera defaultCamera;
// Camera controlled by keyboard and mouse
graphics::FreeCamera freeCamera;
// Current camera
graphics::Camera* currentCamera = &defaultCamera;
// Switch for render camera control panel
bool isUsingCameraPanel = false;
// Is free camera?
bool isUsingFreeCamera = false;
// Mouse is disabled?
bool isMouseBinded = false;
// IK motion
std::unique_ptr<acclaim::Motion> IK;
// Reset motion without reload
std::unique_ptr<acclaim::Motion> IK_backup;
// IK "root" bone
int start_bone = 11;
// IK "touch" bone
int end_bone = 29;
// IK is stable?
bool isStable = true;
}  // namespace

/**
 * @brief When resizing window, we need to update viewport and camera's aspect
 * ratio.
 * @param window current GLFW context
 * @param screenWidth Screen width after resize.
 * @param screenHeight Screen height after resize.
 */
void reshape(GLFWwindow* window, int screenWidth, int screenHeight);

/**
 * @brief Initialize globals and OpenGL
 *
 * @return `GLFWwindow*` The current GLFW context window
 */
GLFWwindow* initialize();
/**
 * @brief Callback function for the debug camera.
 */
void cameraKeyboard(GLFWwindow* window, int key, int scancode, int action, int mode);
/**
 * @brief Reset pointer of std::unique_ptr to call destructor before calling
 * glfwTerminate() to avoid segfault since some destructor contains OpenGL
 * functions.
 *
 * @param window The context to be destroyed
 */
void shutdown();
/**
 * @brief Dear-ImGui main control panel
 *
 */
void mainPanel(kinematics::Ball* ball);
/**
 * @brief Dear-ImGui camera control panel
 *
 * @param window: Current context
 */
void cameraPanel(GLFWwindow* window);
/**
 * @brief Put panels between ImGui::NewFrame() and ImGui::Render() to render
 *
 * @param window: Current context
 */
void renderUI(GLFWwindow* window, kinematics::Ball* ball);

int main() {
    GLFWwindow* window = initialize();
    // No window created
    if (window == nullptr) return 1;
    // Shader programs
    graphics::Program renderProgram;
    graphics::Program skyboxRenderProgram;
    graphics::Program shadowProgram;
    // Texture for shadow mapping
    graphics::ShadowMapTexture shadow(shadowTextureSize);
    // The skybox
    graphics::Box skybox;
    // Ground
    graphics::Plane plane;
    // Marker
    kinematics::Ball ball;
    // Load assets, setup textures
    {
        // Shader
        auto shader_folder = util::PathFinder::find("Shader");
        graphics::Shader shadowVertexShader(shader_folder / "shadow.vert", GL_VERTEX_SHADER);
        graphics::Shader shadowFragmentShader(shader_folder / "shadow.frag", GL_FRAGMENT_SHADER);
        graphics::Shader renderVertexShader(shader_folder / "render.vert", GL_VERTEX_SHADER);
        graphics::Shader renderFragmentShader(shader_folder / "render.frag", GL_FRAGMENT_SHADER);
        graphics::Shader skyboxVertexShader(shader_folder / "skybox.vert", GL_VERTEX_SHADER);
        graphics::Shader skyboxFragmentShader(shader_folder / "skybox.frag", GL_FRAGMENT_SHADER);
        // Setup shaders, these objects can be destroyed after linkShader()
        renderProgram.attachLinkShader(renderVertexShader, renderFragmentShader);
        shadowProgram.attachLinkShader(shadowVertexShader, shadowFragmentShader);
        skyboxRenderProgram.attachLinkShader(skyboxVertexShader, skyboxFragmentShader);
        // Texture
        auto texture_folder = util::PathFinder::find("Texture");
        auto wood = std::make_shared<graphics::Texture>(texture_folder / "wood.png");
        plane.setTexture(std::move(wood));
        auto sky = std::make_shared<graphics::CubeTexture>(std::array<util::fs::path, 6>{
            texture_folder / "skybox0.png", texture_folder / "skybox1.png", texture_folder / "skybox2.png",
            texture_folder / "skybox3.png", texture_folder / "skybox4.png", texture_folder / "skybox5.png"});
        skybox.setTexture(std::move(sky));
        // Acclaim
        auto acclaim_folder = util::PathFinder::find("Acclaim");
        auto skeleton = std::make_unique<acclaim::Skeleton>(acclaim_folder / "skeleton.asf", 0.2);
        IK = std::make_unique<acclaim::Motion>(acclaim_folder / "IK.amc", std::move(skeleton));
        IK_backup = std::make_unique<acclaim::Motion>(*IK);
    }
    // Setup light, uniforms are persisted.
    {
        Eigen::Vector4f lightPosition(11.1f, 24.9f, -14.8f, 0.0f);
        Eigen::Matrix4f lightSpaceMatrix = util::ortho(-30.0f, 30.0f, -30.0f, 30.0f, -75.0f, 75.0f);
        lightSpaceMatrix *= util::lookAt(lightPosition, Eigen::Vector4f::Zero(), Eigen::Vector4f::UnitY());
        // Shader program should be use atleast once before setting up uniforms
        shadowProgram.use();
        shadowProgram.setUniform("lightSpaceMatrix", lightSpaceMatrix);

        renderProgram.use();
        renderProgram.setUniform("lightSpaceMatrix", lightSpaceMatrix);
        renderProgram.setUniform("shadowMap", shadow.getIndex());
        renderProgram.setUniform("lightPos", lightPosition);
    }
    while (!glfwWindowShouldClose(window)) {
        // Moving camera only if debug camera is on.
        if (isUsingFreeCamera) {
            if (isMouseBinded) freeCamera.moveSight(window);
            freeCamera.moveCamera(window);
        }
        currentCamera->update();
        isStable = IK->inverseKinematics(ball.getCurrentPosition(), start_bone, end_bone);
        ball.setModelMatrix();

        // 1. Render shadow to texture
        glViewport(0, 0, shadow.getShadowSize(), shadow.getShadowSize());
        glCullFace(GL_FRONT);
        shadowProgram.use();
        shadow.bindFrameBuffer();
        glClear(GL_DEPTH_BUFFER_BIT);
        plane.render(&shadowProgram);
        ball.render(&shadowProgram);
        IK->render(&shadowProgram);
        shadow.unbindFrameBuffer();
        glCullFace(GL_BACK);
        // 2. Render scene
        glViewport(0, 0, g_ScreenWidth, g_ScreenHeight);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderProgram.use();
        renderProgram.setUniform("viewPos", currentCamera->getPosition());
        renderProgram.setUniform("VP", currentCamera->getViewWithProjectionMatrix());

        plane.render(&renderProgram);
        ball.render(&renderProgram);
        IK->render(&renderProgram);
        // 3. Render the skybox .
        skyboxRenderProgram.use();
        skyboxRenderProgram.setUniform("projection", currentCamera->getProjectionMatrix());
        skyboxRenderProgram.setUniform("view", currentCamera->getViewMatrix());
        skybox.render(&skyboxRenderProgram);
        // 4. Render ImGui UI
        renderUI(window, &ball);
        glFlush();
        glfwSwapBuffers(window);
        // Keyboard and mouse inputs.
        glfwPollEvents();
    }
    shutdown();
    glfwDestroyWindow(window);
    return 0;
}

void reshape(GLFWwindow*, int screenWidth, int screenHeight) {
    g_ScreenWidth = screenWidth;
    g_ScreenHeight = screenHeight;
    glViewport(0, 0, g_ScreenWidth, g_ScreenHeight);
    defaultCamera.setAspectRatio(g_ScreenWidth, g_ScreenHeight);
    freeCamera.setAspectRatio(g_ScreenWidth, g_ScreenHeight);
}

GLFWwindow* initialize() {
    // Initialize GLFW
    if (glfwInit() == GLFW_FALSE) {
        std::cerr << "Failed to initialize GLFW!" << std::endl;
        return nullptr;
    }
    std::atexit(glfwTerminate);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Create GLFW context
    GLFWwindow* window = glfwCreateWindow(g_ScreenWidth, g_ScreenHeight, "Inverse Kinematics", nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create OpenGL 4.1 window!" << std::endl;
        return nullptr;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    // Initialize glad
    if (!gladLoadGL(glfwGetProcAddress)) {
        std::cerr << "Failed to initialize OpenGL context" << std::endl;
        return nullptr;
    }
    // Find assets folder
    if (!util::PathFinder::initialize()) {
        std::cerr << "Cannot find assets!" << std::endl;
        return nullptr;
    }
    // OK, everything works fine
    // ----------------------------------------------------------
    // For high dpi monitors
    glfwGetFramebufferSize(window, &g_ScreenWidth, &g_ScreenHeight);
    GLFWmonitor* moniter = glfwGetPrimaryMonitor();
    const GLFWvidmode* vidMode = glfwGetVideoMode(moniter);
    int maxTextureSize = 1024;
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTextureSize);
    shadowTextureSize = std::min(shadowTextureSize, maxTextureSize);
    // Print some system information
    std::cout << std::left << std::setw(26) << "Current OpenGL renderer"
              << ": " << glGetString(GL_RENDERER) << std::endl;
    std::cout << std::left << std::setw(26) << "Current OpenGL context"
              << ": " << glGetString(GL_VERSION) << std::endl;
    std::cout << std::left << std::setw(26) << "Moniter refresh rate"
              << ": " << vidMode->refreshRate << " Hz" << std::endl;
    std::cout << std::left << std::setw(26) << "Max texture size support"
              << ": " << maxTextureSize << " * " << maxTextureSize << std::endl;
    std::cout << std::left << std::setw(26) << "Shadow texture size"
              << ": " << shadowTextureSize << " * " << shadowTextureSize << std::endl;
    // Setup Opengl
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    // Setup GLFW
    reshape(window, g_ScreenWidth, g_ScreenHeight);
    glfwSetFramebufferSizeCallback(window, reshape);
    // Initialize dear-ImGui
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410 core");
    return window;
}

void cameraKeyboard(GLFWwindow* window, int key, int, int action, int) {
    if (action == GLFW_PRESS && key == GLFW_KEY_F9) {
        if (glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
            // Show the mouse cursor
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            isMouseBinded = false;
        } else {
            // Reset delta x and delta y to avoid view teleporting
            freeCamera.reset();
            // Hide the mouse cursor
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            isMouseBinded = true;
        }
    }
}

void shutdown() {
    IK.reset();
    IK_backup.reset();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void mainPanel(kinematics::Ball* ball) {
    // Main Panel
    ImGui::SetNextWindowSize(ImVec2(300.0f, 210.0f), ImGuiCond_Once);
    ImGui::SetNextWindowCollapsed(0, ImGuiCond_Once);
    ImGui::SetNextWindowPos(ImVec2(60.0f, 525.0f), ImGuiCond_Once);
    ImGui::SetNextWindowBgAlpha(0.2f);
    if (ImGui::Begin("Control")) {
        // Open / close other panels.
        // x ^= 1 === x = !x
        if (ImGui::Button("Camera Panel")) {
            isUsingCameraPanel ^= true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset Skeleton")) {
            IK.reset();
            IK = std::make_unique<acclaim::Motion>(*IK_backup);
        }
        ImGui::SameLine();
        ImGui::Text(isStable ? "Stable" : "Unstable");
        auto&& bpos = ball->getCurrentPosition();
        if (ImGui::InputDouble("target x", &bpos[0], 0.01, 0.1, "%.2lf")) {
            bpos[0] = std::clamp(bpos[0], -10.0, 10.0);
        }
        if (ImGui::InputDouble("target y", &bpos[1], 0.01, 0.1, "%.2lf")) {
            bpos[1] = std::clamp(bpos[1], -10.0, 10.0);
        }
        if (ImGui::InputDouble("target z", &bpos[2], 0.01, 0.1, "%.2lf")) {
            bpos[2] = std::clamp(bpos[2], -10.0, 10.0);
        }
        auto& sk = *IK->getSkeleton();
        ImGui::Text("The last movable bone");
        if (ImGui::InputInt(sk.getBonePointer(start_bone)->name.c_str(), &start_bone, 1, 1)) {
            start_bone = std::clamp(start_bone, 0, sk.getBoneNum() - 1);
        }
        ImGui::Text("The bone that touches the target");
        if (ImGui::InputInt(sk.getBonePointer(end_bone)->name.c_str(), &end_bone, 1, 1)) {
            end_bone = std::clamp(end_bone, 0, sk.getBoneNum() - 1);
        }
    }
    ImGui::End();
}

void cameraPanel(GLFWwindow* window) {
    // Camera control
    // Create imgui window
    ImGui::SetNextWindowSize(ImVec2(500.0f, 210.0f), ImGuiCond_Once);
    ImGui::SetNextWindowCollapsed(0, ImGuiCond_Once);
    ImGui::SetNextWindowPos(ImVec2(370.0f, 525.0f), ImGuiCond_Once);
    ImGui::SetNextWindowBgAlpha(0.2f);
    if (ImGui::Begin("Camera", &isUsingCameraPanel)) {
        if (!isUsingFreeCamera) ImGui::Text("Use this panel to control the camera");
        Eigen::Vector4f cameraPosition = currentCamera->getPosition();
        ImGui::Text("Camera position : (%f, %f, %f)", cameraPosition[0], cameraPosition[1], cameraPosition[2]);
        Eigen::Vector4f cameraCenter = currentCamera->getCenter();
        ImGui::Text("Camera lookat : (%f, %f, %f)", cameraCenter[0], cameraCenter[1], cameraCenter[2]);
        if (!isUsingFreeCamera) {
            ImGui::SliderFloat("Camera rotation angle", defaultCamera.getCameraRotationAnglePointer(), 0.0f, 360.0f);
            ImGui::SliderFloat("Camera rotation radius", defaultCamera.getCameraRotationRadiusPointer(), 0.125f, 50.0f);
            ImGui::SliderFloat("Camera Y Offset", defaultCamera.getCameraYOffsetPointer(), -10.0f, 10.0f);
            if (ImGui::Button("Debug Mode")) {
                isUsingFreeCamera = true;
                currentCamera = &freeCamera;
                glfwSetKeyCallback(window, cameraKeyboard);
            }
        } else {
            ImGui::InputFloat("Mouse sensitivity", freeCamera.getMouseSensitivityPointer(), 0.01f, 0.05f);
            ImGui::InputFloat("Move speed", freeCamera.getMoveSpeedPointer(), 0.01f, 0.05f);
            ImGui::Text("Use W A S D CTRL SPACE to move");
            ImGui::Text("Press F9 to bind / unbind mouse");
            ImGui::Text("Bind mouse to control view");
            if (ImGui::Button("Leave debug Mode")) {
                isUsingFreeCamera = false;
                currentCamera = &defaultCamera;
                glfwSetCursorPosCallback(window, nullptr);
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            }
        }
    }
    ImGui::End();
}

void renderUI(GLFWwindow* window, kinematics::Ball* ball) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    mainPanel(ball);
    if (isUsingCameraPanel) cameraPanel(window);
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
