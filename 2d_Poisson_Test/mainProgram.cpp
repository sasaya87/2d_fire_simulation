// 2d_Poisson_Test.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include "Main.h"
#include "Object.h"
#include "Const.h"


//model行列を定義（中身は単位行列）
mat4 modelMat(1.0);

// View行列を計算
mat4 viewMat = glm::lookAt(
    vec3(dx * l / 2.0, dy * m / 2.0, 3.0), // ワールド空間でのカメラの座標
    vec3(dx * l / 2.0, dy * m / 2.0, 0.0), // 見ている位置の座標
    vec3(0.0, 1.0, 0.0)  // 上方向を示す。(0,1.0,0)に設定するとy軸が上になります
);


mat4 projectionMat = glm::ortho(
    -OrthWidth / 2.0f, //left
    OrthWidth / 2.0f, //right
    -OrthWidth * ((float)Height / (float)Width) / 2.0f, //bottom
    OrthWidth * ((float)Height / (float)Width) / 2.0f, //top
    0.1f,		// 近くのクリッピング平面
    100.0f		// 遠くのクリッピング平面
);

// ModelViewProjection行列を計算
mat4 mvpMat = projectionMat * viewMat * modelMat;



int main()
{
    // GLFW初期化
    if (glfwInit() == GL_FALSE)
    {
        //初期化に失敗した場合
        return -1;
    }

    // ウィンドウ生成
    GLFWwindow* window = glfwCreateWindow(Width, Height, "2D Fire Simulation", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // バージョン3.2指定
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // コンテキストの作成
    glfwMakeContextCurrent(window);

    glfwSwapInterval(1);

    // GLEW初期化
    if (glewInit() != GLEW_OK)
    {
        return -1;
    }

    GLint shader = makeShader("shader.vert", "shader.frag");
    GLuint matrixID = glGetUniformLocation(shader, "MVP");

    //描画に使用するプログラムオブジェクトを指定
    glUseProgram(shader);

    //位置データを作成
    createPositions();

    init();

    std::chrono::system_clock::time_point start, gra1, kawa2, kawa3, kawa0, kawa1, press, velot, gra2;


    //ウィンドウを開いている間繰り返す
    while (glfwWindowShouldClose(window) == GL_FALSE)
    {
        //時間計測をスタート
        start = std::chrono::system_clock::now();
        double time;
        std::cout << "-----TimeStart-----" << std::endl;

        //形状データ作成
        const Object object = createTriangle();

        // デプステストを有効にする
        glEnable(GL_DEPTH_TEST);
        // 前のものよりもカメラに近ければ、フラグメントを受け入れる
        glDepthFunc(GL_LESS);

        // ブレンドの有効化
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // バッファのクリア
        glClearColor(0.2f, 0.2f, 0.2f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // 現在バインドしているシェーダのuniform変数"MVP"に変換行列を送る
        // 4つ目の引数は行列の最初のアドレスを渡しています。
        glUniformMatrix4fv(matrixID, 1, GL_FALSE, &mvpMat[0][0]);

        gra1 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(gra1 - start).count() / 1000.0);
        std::cout << "gra1: " << time << "msec" << std::endl;
        
        //川村スキーム
        omega();
        kawamuraByBiCGSTAB(2);
        kawa2 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(kawa2 - gra1).count() / 1000.0);
        std::cout << "Temperature_advect: " << time << "msec" << std::endl;

        kawamuraByBiCGSTAB(3);
        kawa3 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(kawa3 - kawa2).count() / 1000.0);
        std::cout << "FuelDensity_advect: " << time << "msec" << std::endl;

        kawamuraByBiCGSTAB(0);
        kawa0 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(kawa0 - kawa3).count() / 1000.0);
        std::cout << "U_advect: " << time << "msec" << std::endl;

        kawamuraByBiCGSTAB(1);
        kawa1 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(kawa1 - kawa0).count() / 1000.0);
        std::cout << "V_advect: " << time << "msec" << std::endl;
        
        pressureByCG();
        press = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(press - kawa1).count() / 1000.0);
        std::cout << "Pressure: " << time << "msec" << std::endl;

        velo(true);
        velo(false);
        velot = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(velot - press).count() / 1000.0);
        std::cout << "Velocity_update:" << time << "msec" << std::endl;
        
        //図形の描画
        glBindVertexArray(object.vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, object.count);

        // ダブルバッファのスワップ
        glfwSwapBuffers(window);
        glfwPollEvents();

        gra2 = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(gra2 - press).count() / 1000.0);
        std::cout << "gra2:" << time << "msec" << std::endl;
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(gra2 - start).count() / 1000.0);
        std::cout << "-----TimeEnd_FPS: " << 1000 / time << "fps-----" << std::endl;

    }

    // GLFWの終了処理
    glfwTerminate();

    return 0;
}

// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します
