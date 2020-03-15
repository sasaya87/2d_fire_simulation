#include "Graphics.h"
#include "Const.h"
#include "Object.h"
#include "Physics.h"

//-----------------------------------------------------
//�O���t�B�b�N�X�֘A�̏���
//-----------------------------------------------------


//�ʒu�f�[�^�̔z��
GLfloat position[(2 * l + 4) * m * 2];
//�R���Z�x�f�[�^�̔z��
GLfloat trans[(2 * l + 4) * m];

//�ʒu�f�[�^�̍쐬
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

//�V�F�[�_�[�t�@�C���ǂݍ��݊֐�
int readShaderSource(GLuint shaderObj, std::string fileName)
{
    //�t�@�C���̓ǂݍ���
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

    // �V�F�[�_�̃\�[�X�v���O�������V�F�[�_�I�u�W�F�N�g�֓ǂݍ���
    const GLchar* sourcePtr = (const GLchar*)source.c_str();
    GLint length = source.length();
    glShaderSource(shaderObj, 1, &sourcePtr, &length);

    return 0;
}

GLint makeShader(std::string vertexFileName, std::string fragmentFileName)
{
    // �V�F�[�_�[�I�u�W�F�N�g�쐬
    GLuint vertShaderObj = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragShaderObj = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint shader;

    // �V�F�[�_�[�R���p�C���ƃ����N�̌��ʗp�ϐ�
    GLint compiled, linked;

    /* �V�F�[�_�[�̃\�[�X�v���O�����̓ǂݍ��� */
    if (readShaderSource(vertShaderObj, vertexFileName)) return -1;
    if (readShaderSource(fragShaderObj, fragmentFileName)) return -1;

    /* �o�[�e�b�N�X�V�F�[�_�[�̃\�[�X�v���O�����̃R���p�C�� */
    glCompileShader(vertShaderObj);
    glGetShaderiv(vertShaderObj, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE)
    {
        fprintf(stderr, "Compile error in vertex shader.\n");
        return -1;
    }

    /* �t���O�����g�V�F�[�_�[�̃\�[�X�v���O�����̃R���p�C�� */
    glCompileShader(fragShaderObj);
    glGetShaderiv(fragShaderObj, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE)
    {
        fprintf(stderr, "Compile error in fragment shader.\n");
        return -1;
    }

    /* �v���O�����I�u�W�F�N�g�̍쐬 */
    shader = glCreateProgram();

    /* �V�F�[�_�[�I�u�W�F�N�g�̃V�F�[�_�[�v���O�����ւ̓o�^ */
    glAttachShader(shader, vertShaderObj);
    glAttachShader(shader, fragShaderObj);

    /* �V�F�[�_�[�I�u�W�F�N�g�̍폜 */
    glDeleteShader(vertShaderObj);
    glDeleteShader(fragShaderObj);

    /* �V�F�[�_�[�v���O�����̃����N */
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

//���_�z��I�u�W�F�N�g�̍쐬
//***vertices:���_�̐�
//***position:���_�̈ʒu���i�[�����z��
GLuint createObject(GLuint vertices, const GLfloat* position, const GLfloat* trans)
{
    //���_�z��I�u�W�F�N�g
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    //���_�o�b�t�@�I�u�W�F�N�g(�ʒu)
    GLuint vboPos;
    glGenBuffers(1, &vboPos);
    glBindBuffer(GL_ARRAY_BUFFER, vboPos);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * 2 * vertices, position, GL_STATIC_DRAW);

    //��������Ă��钸�_�o�b�t�@�I�u�W�F�N�g��attribute�ϐ�����Q�Ƃł���悤�ɂ���
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    //���_�o�b�t�@�I�u�W�F�N�g(�����x)
    GLuint vboClr;
    glGenBuffers(1, &vboClr);
    glBindBuffer(GL_ARRAY_BUFFER, vboClr);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * 1 * vertices, trans, GL_STATIC_DRAW);

    //��������Ă��钸�_�o�b�t�@�I�u�W�F�N�g��attribute�ϐ�����Q�Ƃł���悤�ɂ���
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);

    //VBO��VAO�̌������������
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return vao;
}


//�O�p�`�̃f�[�^���쐬����
Object createTriangle() {
    

    //�R���Z�x�f�[�^
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
    

    


    //���_�̐�
    static const int vertices = sizeof(position) / (2 * sizeof(position[0]));

    //���_�z��I�u�W�F�N�g�̍쐬
    Object object;
    object.vao = createObject(vertices, position, trans);
    object.count = vertices;

    return object;
}
