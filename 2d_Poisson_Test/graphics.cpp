#include "Graphics.h"
#include "Const.h"
#include "Object.h"
#include "Physics.h"

//-----------------------------------------------------
//グラフィックス関連の処理
//-----------------------------------------------------


//位置データの配列
GLfloat position[(2 * l + 4) * m * 2];
//燃料濃度データの配列
GLfloat trans[(2 * l + 4) * m];

//位置データの作成
void createPositions() {

    for (int j = 0; j < m; j++) {
        for (int i = 0; i <= l; i++) {
            if (i == 0) {
                position[((2 * l + 4) * j + 2 * i) * 2] = (GLfloat)dx * i;
                position[((2 * l + 4) * j + 2 * i) * 2 + 1] = (GLfloat)dy * j;
            }
            position[((2 * l + 4) * j + 2 * i + 1) * 2] = (GLfloat)dx * i;
            position[((2 * l + 4) * j + 2 * i + 1) * 2 + 1] = (GLfloat)dy * j;
            position[((2 * l + 4) * j + 2 * i + 2) * 2] = (GLfloat)dx * i;
            position[((2 * l + 4) * j + 2 * i + 2) * 2 + 1] = (GLfloat)dy * (j + 1);
            if (i == l) {
                position[((2 * l + 4) * j + 2 * i + 3) * 2] = (GLfloat)dx * i;
                position[((2 * l + 4) * j + 2 * i + 3) * 2 + 1] = (GLfloat)dy * (j + 1);
            }
        }
    }
}

//シェーダーファイル読み込み関数
int readShaderSource(GLuint shaderObj, std::string fileName)
{
    //ファイルの読み込み
    std::ifstream ifs(fileName);
    if (!ifs)
    {
        std::cout << "error" << std::endl;
        return -1;
    }

    std::string source;
    std::string line;
    while (getline(ifs, line))
    {
        source += line + "\n";
    }

    // シェーダのソースプログラムをシェーダオブジェクトへ読み込む
    const GLchar* sourcePtr = (const GLchar*)source.c_str();
    GLint length = source.length();
    glShaderSource(shaderObj, 1, &sourcePtr, &length);

    return 0;
}

GLint makeShader(std::string vertexFileName, std::string fragmentFileName)
{
    // シェーダーオブジェクト作成
    GLuint vertShaderObj = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragShaderObj = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint shader;

    // シェーダーコンパイルとリンクの結果用変数
    GLint compiled, linked;

    /* シェーダーのソースプログラムの読み込み */
    if (readShaderSource(vertShaderObj, vertexFileName)) return -1;
    if (readShaderSource(fragShaderObj, fragmentFileName)) return -1;

    /* バーテックスシェーダーのソースプログラムのコンパイル */
    glCompileShader(vertShaderObj);
    glGetShaderiv(vertShaderObj, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE)
    {
        fprintf(stderr, "Compile error in vertex shader.\n");
        return -1;
    }

    /* フラグメントシェーダーのソースプログラムのコンパイル */
    glCompileShader(fragShaderObj);
    glGetShaderiv(fragShaderObj, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE)
    {
        fprintf(stderr, "Compile error in fragment shader.\n");
        return -1;
    }

    /* プログラムオブジェクトの作成 */
    shader = glCreateProgram();

    /* シェーダーオブジェクトのシェーダープログラムへの登録 */
    glAttachShader(shader, vertShaderObj);
    glAttachShader(shader, fragShaderObj);

    /* シェーダーオブジェクトの削除 */
    glDeleteShader(vertShaderObj);
    glDeleteShader(fragShaderObj);

    /* シェーダープログラムのリンク */
    glBindAttribLocation(shader, 0, "position");
    glBindAttribLocation(shader, 1, "trans");
    glBindFragDataLocation(shader, 0, "fragColor");
    glLinkProgram(shader);
    glGetProgramiv(shader, GL_LINK_STATUS, &linked);
    if (linked == GL_FALSE)
    {
        fprintf(stderr, "Link error.\n");
        return -1;
    }

    return shader;
}

//頂点配列オブジェクトの作成
//***vertices:頂点の数
//***position:頂点の位置を格納した配列
GLuint createObject(GLuint vertices, const GLfloat* position, const GLfloat* trans)
{
    //頂点配列オブジェクト
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    //頂点バッファオブジェクト(位置)
    GLuint vboPos;
    glGenBuffers(1, &vboPos);
    glBindBuffer(GL_ARRAY_BUFFER, vboPos);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * 2 * vertices, position, GL_STATIC_DRAW);

    //結合されている頂点バッファオブジェクトをattribute変数から参照できるようにする
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    //頂点バッファオブジェクト(透明度)
    GLuint vboClr;
    glGenBuffers(1, &vboClr);
    glBindBuffer(GL_ARRAY_BUFFER, vboClr);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * 1 * vertices, trans, GL_STATIC_DRAW);

    //結合されている頂点バッファオブジェクトをattribute変数から参照できるようにする
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);

    //VBOとVAOの結合を解放する
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return vao;
}


//三角形のデータを作成する
Object createTriangle() {
    

    //燃料濃度データ
    for (int j = 0; j < m; j++) {
        for (int i = 0; i <= l; i++) {
            if (i == 0) {
                trans[(2 * l + 4) * j + 2 * i] = rho[i][j] / rhob;
            }
            trans[(2 * l + 4) * j + 2 * i + 1] = rho[i][j] / rhob;
            trans[(2 * l + 4) * j + 2 * i + 2] = rho[i][j + 1] / rhob;
            if (i == l) {
                trans[(2 * l + 4) * j + 2 * i + 3] = rho[i][j + 1] / rhob;
            }
        }
    }
    

    


    //頂点の数
    static const int vertices = sizeof(position) / (2 * sizeof(position[0]));

    //頂点配列オブジェクトの作成
    Object object;
    object.vao = createObject(vertices, position, trans);
    object.count = vertices;

    return object;
}
