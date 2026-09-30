#include <raylib.h>

int main()
{
    const int screenWidth = 1024;
    const int screenHeight = 768;
    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera mode");
    Camera3D camera = {};
    camera.position = (Vector3){ 0.0f, 10.0f, 10.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    Vector3 cylinderPosition = { 0.0f, 0.0f, 0.0f };
    Model cylinder = LoadModelFromMesh(GenMeshCylinder(1.0f, 2.0f, 32));
    float rotation = 0.0f;
    SetTargetFPS(60);
    while(!WindowShouldClose())
    {
        rotation += 1.0f;
        BeginDrawing();
        ClearBackground(BLACK);
        BeginMode3D(camera);
        DrawModelWiresEx(cylinder, cylinderPosition, (Vector3){ 0.0f, 1.0f, 0.0f }, rotation, (Vector3){ 2.0f, 2.0f, 2.0f }, GRAY);
        EndMode3D();
        EndDrawing();
    }
    UnloadModel(cylinder);
    CloseWindow();
    return 0;
}
