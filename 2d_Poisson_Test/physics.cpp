#include <algorithm>
#include <iostream>
#include <omp.h> 
#include "Physics.h"

double rho[l + 3][m + 3]; //�R���̖��x

double tempera[l + 3][m + 3]; //���x

double u3[l + 3][m + 3]; //�ڗ��g�U�v�Z���u�̒��ԑ��x
double v3[l + 3][m + 3]; //�ڗ��g�U�v�Z���v�̒��ԑ��x

double u[l + 3][m + 3]; //x�������xu
double v[l + 3][m + 3]; //y�������xv

double ome[l + 3][m + 3]; //�Q�x

double p[l + 3][m + 3]; //����
double px[l + 3][m + 3]; //x�����̈��͌��z
double py[l + 3][m + 3]; //y�����̈��͌��z

//���͂̃|�A�\���������̌W��
double cp1 = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
double cp2 = 1.0 / (dx * dx);
double cp3 = 1.0 / (dx * dx);
double cp4 = 1.0 / (dy * dy);
double cp5 = 1.0 / (dy * dy);

//�ڗ��g�U�������̌W��
double cx1[l + 3][m + 3];
double cx2[l + 3][m + 3];
double cx3[l + 3][m + 3];
double cx4[l + 3][m + 3];
double cx5[l + 3][m + 3];
double cx6[l + 3][m + 3];
double cx7[l + 3][m + 3];
double cx8[l + 3][m + 3];
double cx9[l + 3][m + 3];


//-----------------------------------------------------------
//����������
//-----------------------------------------------------------
void init() {
    
    int i, j;

    //���x�̏�����
    #pragma omp parallel
    {
        #pragma omp for private(j)
        for (i = 2; i < l + 1; i++) {
            for (j = 2; j < m / 2 + 1; j++) {
                tempera[i][j] = Tmin;
            }
        }
        #pragma omp for
        for (i = 0; i < l + 3; i++) {
            tempera[i][0] = Tmin;
            tempera[i][1] = Tmin;
            tempera[i][m + 1] = Tmin;
            tempera[i][m + 2] = Tmin;
        }
        #pragma omp for
        for (j = 0; j < m + 3; j++) {
            tempera[0][j] = Tmin;
            tempera[1][j] = Tmin;
            tempera[l + 1][j] = Tmin;
            tempera[l + 2][j] = Tmin;
        }
        #pragma omp for
        for (j = 9 * m / 20; j < 11 * m / 20; j++) {
            tempera[1][j] = Tmax;
        }

        //���x�A���́A���͂̏�����
        #pragma omp for private(j)
        for (i = 0; i < l + 3; i++) {
            for (j = 0; j < m + 3; j++) {
                u[i][j] = 0.0;
                v[i][j] = 0.0;
                u3[i][j] = 0.0;
                v3[i][j] = 0.0;
                p[i][j] = pb;
                px[i][j] = pb;
                py[i][j] = pb;
                rho[i][j] = 0.0;
                ome[i][j] = 0.0;
                cx1[i][j] = 0.0;
                cx2[i][j] = 0.0;
                cx3[i][j] = 0.0;
                cx4[i][j] = 0.0;
                cx5[i][j] = 0.0;
                cx6[i][j] = 0.0;
                cx7[i][j] = 0.0;
                cx8[i][j] = 0.0;
                cx9[i][j] = 0.0;
            }
        }

        #pragma omp for
        for (j = 9 * m / 20; j < 11 * m / 20; j++) {
            u[1][j] = ub; //�����o������u��ub
            v[1][j] = vb; //�����o������v��vb
            rho[1][j] = rhob; //�����o������rho��rhob
        }
    }
    std::cout << "init" << std::endl;
}

//-----------------------------------------------------------
//�Q�x�����@
//-----------------------------------------------------------
//�Q�x(��z����)���v�Z(2�����Ȃ̂�z�����̂�)
void omega() {
    int i, j;
    #pragma omp parallel for private(j)
    for ( i = 1; i < l + 2; i++) {
        for ( j = 1; j < m + 2; j++) {
            ome[i][j] = ((v[i + 1][j] - v[i - 1][j]) / dx - (u[i][j + 1] - u[i][j - 1]) / dy) / 2.0;
        }
    }
}
double vortCon(int i, int j, bool isU) {
    double etax = (std::abs(ome[i + 1][j]) - std::abs(ome[i - 1][j])) / (2.0 * dx); //�ł�x����
    double etay = (std::abs(ome[i][j + 1]) - std::abs(ome[i][j - 1])) / (2.0 * dy); //�ł�y����
    double eta = pow(etax * etax + etay * etay, 0.5) + 1.0e-5; //�ł̑傫��
    double nx = etax / eta;
    double ny = etay / eta;
    if (isU) {
        return come * dx * ny * ome[i][j];
    }
    else {
        return -come * dy * nx * ome[i][j];
    }
}

//-----------------------------------------------------------
//�ڗ��g�U�{�O�͍��̌v�Z(�쑺�X�L�[���ABiCGSTAB�@)
//-----------------------------------------------------------
void kawamuraByBiCGSTAB(int mode) {
    //mode=0�̎���u�A1�̎���v�A2�̎���T�A3�̎���rho�̈ڗ��g�U���������v�Z
    double(*ans)[m + 3]; //u,v,tempera,rho�̂ǂ��炩
    double(*ans3)[m + 3]; //u3,v3,tempera,rho�̂ǂ��炩
    double tempc; //a��1/re
    double tempp; //�����o�����ȊO��T�A���邢��rho
    double tempb; //�����o������T�A���邢��rho

    bool continueFlag = true;
    int cnt = 0;

    if (mode == 0) {
        ans = u;
        ans3 = u3;
        tempc = 1 / re;
        tempp = 0.0;
        tempb = ub;
    }
    else if (mode == 1) {
        ans = v;
        ans3 = v3;
        tempc = 1 / re;
        tempp = 0.0;
        tempb = vb;
    }
    else if (mode == 2) {
        ans = tempera;
        ans3 = tempera;
        tempc = a;
        tempp = Tmin;
        tempb = Tmax;
    }
    else {
        ans = rho;
        ans3 = rho;
        tempc = arho;
        tempp = 0.0;
        tempb = rhob;
    }

    double resi[l + 3][m + 3]; //�_(i,j)�ɂ�����c���x�N�g��r
    double resia[l + 3][m + 3]; //�_(i,j)�ɂ�����V���h�E�c���x�N�g��r*
    double pp[l + 3][m + 3]; //�_(i,j)�ɂ�����C�������x�N�g��p
    double s[l + 3][m + 3]; 
    double y[l + 3][m + 3]; //App
    double z[l + 3][m + 3]; //As
    double alpha = 0.0; //���̏C���W����
    double beta = 0.0; //�����x�N�g��p�̏C���W����
    double omg = 0.0; //�C���W����
    double r2 = 0.0; //�c��resi���m�̓���
    double rr = 0.0; //�V���h�E�c��resia�Ə����c��r0�̓���
    double nrr = 0.0; //���̃X�e�b�v�̎c��resi�ƃV���h�E�c��resia�̓���
    double r0y = 0.0; //r0*��y=App�̓���
    double zs = 0.0; //z��s�̓���
    double zz = 0.0; //z���m�̓���

    int i, j;

    //�W���̌v�Z(�ڗ��������̌v�Z�̎��̂�)
    #pragma omp parallel
    {
        #pragma omp for private(j)
        for (i = 0; i < l + 3; i++) {
            for (j = 0; j < m + 3; j++) {
                cx1[i][j] = 1.0 / dt + std::abs(u[i][j]) * 3.0 / (2.0 * dx) + std::abs(v[i][j]) * 3.0 / (2.0 * dy)
                    + tempc * 2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)); //u3[i][j]�̌W��
                cx2[i][j] = -u[i][j] * 2.0 / (3.0 * dx) - std::abs(u[i][j]) / dx - tempc / (dx * dx); //u3[i-1][j]�̌W��
                cx3[i][j] = u[i][j] * 2.0 / (3.0 * dx) - std::abs(u[i][j]) / dx - tempc / (dx * dx); //u3[i+1][j]�̌W��
                cx4[i][j] = -v[i][j] * 2.0 / (3.0 * dy) - std::abs(v[i][j]) / dy - tempc / (dy * dy); //u3[i][j-1]�̌W��
                cx5[i][j] = v[i][j] * 2.0 / (3.0 * dy) - std::abs(v[i][j]) / dy - tempc / (dy * dy); //u3[i][j+1]�̌W��
                cx6[i][j] = u[i][j] / (12.0 * dx) + std::abs(u[i][j]) / (4.0 * dx); //u3[i-2][j]�̌W��
                cx7[i][j] = -u[i][j] / (12.0 * dx) + std::abs(u[i][j]) / (4.0 * dx); //u3[i+2][j]�̌W��
                cx8[i][j] = v[i][j] / (12.0 * dy) + std::abs(v[i][j]) / (4.0 * dy); //u3[i][j-2]�̌W��
                cx9[i][j] = -v[i][j] / (12.0 * dy) + std::abs(v[i][j]) / (4.0 * dy); //u3[i][j+2]�̌W��
            }
        }


        //BiCG�@�̏����ݒ�
        #pragma omp for private(j) reduction(+:rr)
        for (i = 2; i < l + 1; i++) {
            for (j = 2; j < m + 1; j++) {
                double f = 0.0; //�O�́A�������̓\�[�X��
                if (mode == 1) {
                    //v�̏ꍇ�̂ݕ��͂�������(y�����̕���)
                    //f = beta * (tempera[i][j] - Tmin); //���͂ƉQ
                    f = beta * (tempera[i][j] - Tmin) + vortCon(i, j, false); //���͂ƉQ
                }
                else if (mode == 3) {
                    //�R�����x�͎��Ԍo�߂ɏ]���Ďw���֐��I�Ɍ���
                    f = -crho * rho[i][j];
                }
                else if (mode == 2) {
                    //�R���̖��x�ɔ�Ⴕ�ĔM�𔭐�
                    f = rhot * rho[i][j];
                }
                else if (mode == 0) {
                    f = vortCon(i, j, true);
                }

                double xcon = ans[i][j] / dt + f; //cx1�Ŋ���O�̉E��

                resi[i][j] = xcon - (cx1[i][j] * ans3[i][j] + cx2[i][j] * ans3[i - 1][j] + cx3[i][j] * ans3[i + 1][j]
                    + cx4[i][j] * ans3[i][j - 1] + cx5[i][j] * ans3[i][j + 1]
                    + cx6[i][j] * ans3[i - 2][j] + cx7[i][j] * ans3[i + 2][j]
                    + cx8[i][j] * ans3[i][j - 2] + cx9[i][j] * ans3[i][j + 2]);

                resia[i][j] = resi[i][j];
                pp[i][j] = resi[i][j];
                rr += resia[i][j] * resi[i][j];
            }
        }

        #pragma omp for
        for (i = 0; i < l + 3; i++) {
            pp[i][0] = 0;
            pp[i][1] = 0;
            pp[i][m + 1] = 0;
            pp[i][m + 2] = 0;
            resi[i][0] = 0;
            resi[i][1] = 0;
            resi[i][m + 1] = 0;
            resi[i][m + 2] = 0;
            resia[i][0] = 0;
            resia[i][1] = 0;
            resia[i][m + 1] = 0;
            resia[i][m + 2] = 0;
            s[i][0] = 0;
            s[i][1] = 0;
            s[i][m + 1] = 0;
            s[i][m + 2] = 0;
        }

        #pragma omp for
        for (j = 0; j < m + 3; j++) {
            resi[0][j] = 0;
            resi[1][j] = 0;
            resi[l + 1][j] = 0;
            resi[l + 2][j] = 0;
            resia[0][j] = 0;
            resia[1][j] = 0;
            resia[l + 1][j] = 0;
            resia[l + 2][j] = 0;
            pp[0][j] = 0;
            pp[1][j] = 0;
            pp[l + 1][j] = 0;
            pp[l + 2][j] = 0;
            s[0][j] = 0;
            s[1][j] = 0;
            s[l + 1][j] = 0;
            s[l + 2][j] = 0;
        }
    }

    while (continueFlag) {
        //y=App���v�Z
        #pragma omp parallel for private(j) reduction(+:r0y)
        for ( i = 2; i < l + 1; i++) {
            for ( j = 2; j < m + 1; j++) {
                if (mode == 4) {
                    y[i][j] = cp3 * pp[i + 1][j] + cp2 * pp[i - 1][j] + cp5 * pp[i][j + 1] + cp4 * pp[i][j - 1] + cp1 * pp[i][j];
                }
                else {
                    y[i][j] = cx1[i][j] * pp[i][j] + cx2[i][j] * pp[i - 1][j] + cx3[i][j] * pp[i + 1][j]
                        + cx4[i][j] * pp[i][j - 1] + cx5[i][j] * pp[i][j + 1]
                        + cx6[i][j] * pp[i - 2][j] + cx7[i][j] * pp[i + 2][j]
                        + cx8[i][j] * pp[i][j - 2] + cx9[i][j] * pp[i][j + 2];
                }
                r0y += resia[i][j] * y[i][j];
            }
        }

        //�C���W�������v�Z
        alpha = rr / r0y;

        #pragma omp parallel
        {
            //s���v�Z
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    s[i][j] = resi[i][j] - alpha * y[i][j];
                }
            }

            //z=As���v�Z
            #pragma omp for private(j) reduction(+:zs) reduction(+:zz)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    z[i][j] = cx1[i][j] * s[i][j] + cx2[i][j] * s[i - 1][j] + cx3[i][j] * s[i + 1][j]
                        + cx4[i][j] * s[i][j - 1] + cx5[i][j] * s[i][j + 1]
                        + cx6[i][j] * s[i - 2][j] + cx7[i][j] * s[i + 2][j]
                        + cx8[i][j] * s[i][j - 2] + cx9[i][j] * s[i][j + 2];
                    zs += z[i][j] * s[i][j];
                    zz += z[i][j] * z[i][j];
                }
            }
        }

        //�C���W���ւ��v�Z
        omg = zs / zz;

        #pragma omp parallel
        {
            //���̃X�e�b�v�̋ߎ��l���v�Z
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    ans3[i][j] += alpha * pp[i][j] + omg * s[i][j];
                }
            }

            //���̃X�e�b�v�̋ߎ��ɑ΂���c�����v�Z
            #pragma omp for private(j) reduction(+:r2) reduction(+:nrr)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    resi[i][j] = s[i][j] - omg * z[i][j];
                    r2 += resi[i][j] * resi[i][j];
                    nrr += resia[i][j] * resi[i][j];
                }
            }
        }

        //�c�����������Ȃ邩���[�v�������l�ȏ�ɂȂ�ΏI��
        if (pow(r2, 0.5) < 1.0e-5 || cnt > l* m) {
            continueFlag = false;
        }

        //beta���v�Z
        beta = (alpha /omg) * (nrr / rr);

        //���̃X�e�b�v�̕����x�N�g��pp���v�Z
        #pragma omp parallel for private(j)
        for ( i = 2; i < l + 1; i++) {
            for ( j = 2; j < m + 1; j++) {
                pp[i][j] = resi[i][j] + beta * (pp[i][j]-omg * y[i][j]);
            }
        }

        rr = nrr;
        nrr = 0.0;
        r2 = 0.0;
        r0y = 0.0;
        zs = 0.0;
        zz = 0.0;

        cnt++;
    }

    //���E������ݒ�
    #pragma omp parallel
    {
        #pragma omp for
        for (i = 0; i < l + 3; i++) {
            ans3[i][0] = ans3[i][2]; //���z�Z��
            ans3[i][1] = tempp; //���̖ʂ̉��xTmin
            ans3[i][m + 1] = tempp; //��̖ʂ̉��xTmin
            ans3[i][m + 2] = ans3[i][m]; //���z�Z��
        }
        #pragma omp for
        for (j = 0; j < m + 3; j++) {
            ans3[0][j] = ans3[2][j]; //���z�Z��
            ans3[1][j] = tempp; //���̖ʂ̉��xTmin
            ans3[l + 1][j] = tempp; //���̖ʂ̉��xTmin
            ans3[l + 2][j] = ans3[l][j]; //���z�Z��
        }
        #pragma omp for
        for (j = 9 * m / 20; j < 11 * m / 20; j++) {
            ans3[1][j] = tempb; //�����o�����̉��xTmax
        }
    }
}

//-----------------------------------------------------
//���͍��̌v�Z(CG�@)
//-----------------------------------------------------
void pressureByCG() {
    bool continueFlag = true;
    int cnt = 0;

    double h; //���͂̃|�A�\���������̉E��

    double resi[l + 3][m + 3]; //�_(i,j)�ɂ�����c��
    double pp[l + 3][m + 3]; //�_(i,j)�ɂ�����C�������x�N�g��
    double y[l + 3][m + 3]; //App
    double alpha = 0.0; //���̏C���W����
    double beta = 0.0; //�����x�N�g��p�̏C���W����
    double rr = 0.0; //�c��resi�̓���
    double nrr = 0.0; //���̃X�e�b�v�̎c��resi�̓���
    double ppy = 0.0; //pp��y�̓���

    int i, j;

    //CG�@�̏����ݒ�
    #pragma omp parallel
    {
        #pragma omp for private(j) reduction(+:rr)
        for (i = 2; i < l + 1; i++) {
            for (j = 2; j < m + 1; j++) {
                h = ((u3[i + 1][j] - u3[i - 1][j]) / dx + (v3[i][j + 1] - v3[i][j - 1]) / dy) / (2 * dt); //���͂̃|�A�\���������̉E��
                resi[i][j] = h - (cp3 * p[i + 1][j] + cp2 * p[i - 1][j] + cp5 * p[i][j + 1] + cp4 * p[i][j - 1] + cp1 * p[i][j]);
                pp[i][j] = resi[i][j];
                rr += resi[i][j] * resi[i][j];
            }
        }
        #pragma omp for
        for (i = 0; i < l + 3; i++) {
            pp[i][0] = 0;
            pp[i][1] = 0;
            pp[i][m + 1] = 0;
            pp[i][m + 2] = 0;
            resi[i][0] = 0;
            resi[i][1] = 0;
            resi[i][m + 1] = 0;
            resi[i][m + 2] = 0;
        }
        #pragma omp for
        for (j = 0; j < m + 3; j++) {
            resi[0][j] = 0;
            resi[1][j] = 0;
            resi[l + 1][j] = 0;
            resi[l + 2][j] = 0;
            pp[0][j] = 0;
            pp[1][j] = 0;
            pp[l + 1][j] = 0;
            pp[l + 2][j] = 0;
        }
    }

    while (continueFlag) {

        //y���v�Z
        #pragma omp parallel for private(j) reduction(+:ppy)
        for ( i = 2; i < l + 1; i++) {
            for ( j = 2; j < m + 1; j++) {
                y[i][j] = cp1 * pp[i][j] + cp2 * pp[i - 1][j] + cp3 * pp[i + 1][j] + cp4 * pp[i][j - 1] + cp5 * pp[i][j + 1];
                ppy += pp[i][j] * y[i][j];
            }
        }

        //�C���W�������v�Z
        alpha = rr / ppy;

        #pragma omp parallel
        {
            //���̃X�e�b�v�̋ߎ��l���v�Z
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    p[i][j] += alpha * pp[i][j];
                }
            }

            //���̃X�e�b�v�̋ߎ��ɑ΂���c�����v�Z
            #pragma omp for private(j) reduction(+:nrr)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    resi[i][j] -= alpha * y[i][j];
                    nrr += resi[i][j] * resi[i][j];
                }
            }
        }

        //�c�����������Ȃ邩���[�v�������l�ȏ�ɂȂ�ΏI��
        if (pow(nrr, 0.5) < 1.0e-5 || cnt > l* m) {
            continueFlag = false;
        }

        //beta���v�Z
        beta = nrr / rr;

        //���̃X�e�b�v�̕����x�N�g��pp���v�Z
        #pragma omp parallel for private(j)
        for ( i = 2; i < l + 1; i++) {
            for ( j = 2; j < m + 1; j++) {
                pp[i][j] = resi[i][j] + beta * pp[i][j];
            }
        }

        rr = nrr;
        ppy = 0.0;
        nrr = 0.0;

        cnt++;
    }

    //���͌��z�̌v�Z
    #pragma omp parallel for private(j)
    for ( i = 1; i < l + 2; i++) {
        for ( j = 1; j < m + 2; j++) {
            px[i][j] = (p[i + 1][j] - p[i - 1][j]) / (2.0 * dx);
            py[i][j] = (p[i][j + 1] - p[i][j - 1]) / (2.0 * dy);
        }
    }
}

//-----------------------------------------------------
//���x�̍X�V
//-----------------------------------------------------
void velo(bool isU) {
    double(*ans)[m + 3]; //u,v�̂ǂ��炩
    double(*ans3)[m + 3]; //u3,v3�̂ǂ��炩
    double(*ansp)[m + 3]; //px,py�̂ǂ��炩
    double tempb;

    if (isU) {
        ans = u;
        ans3 = u3;
        ansp = px;
        tempb = ub;
    }
    else {
        ans = v;
        ans3 = v3;
        ansp = py;
        tempb = vb;
    }

    int i, j;
    #pragma omp parallel 
    {
        #pragma omp for private(j)
        for (i = 1; i < l + 2; i++) {
            for (j = 1; j < m + 2; j++) {
                ans[i][j] = ans3[i][j] - dt * ansp[i][j];
            }
        }


        //���x�̋��E����
        #pragma omp for
        for (int i = 0; i < l + 3; i++) {
            ans[i][0] = ans[i][2]; //���z�Z��
            ans[i][1] = 0.0; //���E�ʂ̑��x0
            ans[i][m + 1] = 0.0; //���E�ʂ̑��x0
            ans[i][m + 2] = ans[i][m]; //���z�Z��
        }
        #pragma omp for
        for (int j = 0; j < m + 3; j++) {
            ans[0][j] = ans[2][j]; //���z�Z��
            ans[1][j] = 0.0; //���E�ʂ̑��x0
            ans[l + 1][j] = 0.0; //���E�ʂ̑��x0
            ans[l + 2][j] = ans[l][j]; //���z�Z��
        }
        #pragma omp for
        for (int j = 9 * m / 20; j < 11 * m / 20; j++) {
            ans[1][j] = tempb; //�����o������u��ub�Av��vb
        }
    }
}
