#include <algorithm>
#include <iostream>
#include <omp.h> 
#include "Physics.h"

double rho[l + 3][m + 3]; //燃料の密度

double tempera[l + 3][m + 3]; //温度

double u3[l + 3][m + 3]; //移流拡散計算後のuの中間速度
double v3[l + 3][m + 3]; //移流拡散計算後のvの中間速度

double u[l + 3][m + 3]; //x方向速度u
double v[l + 3][m + 3]; //y方向速度v

double ome[l + 3][m + 3]; //渦度

double p[l + 3][m + 3]; //圧力
double px[l + 3][m + 3]; //x方向の圧力勾配
double py[l + 3][m + 3]; //y方向の圧力勾配

//圧力のポアソン方程式の係数
double cp1 = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
double cp2 = 1.0 / (dx * dx);
double cp3 = 1.0 / (dx * dx);
double cp4 = 1.0 / (dy * dy);
double cp5 = 1.0 / (dy * dy);

//移流拡散方程式の係数
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
//初期化処理
//-----------------------------------------------------------
void init() {
    
    int i, j;

    //温度の初期化
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

        //速度、圧力、浮力の初期化
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
            u[1][j] = ub; //吹き出し口のuはub
            v[1][j] = vb; //吹き出し口のvはvb
            rho[1][j] = rhob; //吹き出し口のrhoはrhob
        }
    }
    std::cout << "init" << std::endl;
}

//-----------------------------------------------------------
//渦度強制法
//-----------------------------------------------------------
//渦度(のz成分)を計算(2次元なのでz成分のみ)
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
    double etax = (std::abs(ome[i + 1][j]) - std::abs(ome[i - 1][j])) / (2.0 * dx); //ηのx成分
    double etay = (std::abs(ome[i][j + 1]) - std::abs(ome[i][j - 1])) / (2.0 * dy); //ηのy成分
    double eta = pow(etax * etax + etay * etay, 0.5) + 1.0e-5; //ηの大きさ
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
//移流拡散＋外力項の計算(川村スキーム、BiCGSTAB法)
//-----------------------------------------------------------
void kawamuraByBiCGSTAB(int mode) {
    //mode=0の時はu、1の時はv、2の時はT、3の時はrhoの移流拡散方程式を計算
    double(*ans)[m + 3]; //u,v,tempera,rhoのどちらか
    double(*ans3)[m + 3]; //u3,v3,tempera,rhoのどちらか
    double tempc; //aか1/re
    double tempp; //吹き出し口以外のT、あるいはrho
    double tempb; //吹き出し口のT、あるいはrho

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

    double resi[l + 3][m + 3]; //点(i,j)における残差ベクトルr
    double resia[l + 3][m + 3]; //点(i,j)におけるシャドウ残差ベクトルr*
    double pp[l + 3][m + 3]; //点(i,j)における修正方向ベクトルp
    double s[l + 3][m + 3]; 
    double y[l + 3][m + 3]; //App
    double z[l + 3][m + 3]; //As
    double alpha = 0.0; //解の修正係数α
    double beta = 0.0; //方向ベクトルpの修正係数β
    double omg = 0.0; //修正係数ω
    double r2 = 0.0; //残差resi同士の内積
    double rr = 0.0; //シャドウ残差resiaと初期残差r0の内積
    double nrr = 0.0; //次のステップの残差resiとシャドウ残差resiaの内積
    double r0y = 0.0; //r0*とy=Appの内積
    double zs = 0.0; //zとsの内積
    double zz = 0.0; //z同士の内積

    int i, j;

    //係数の計算(移流方程式の計算の時のみ)
    #pragma omp parallel
    {
        #pragma omp for private(j)
        for (i = 0; i < l + 3; i++) {
            for (j = 0; j < m + 3; j++) {
                cx1[i][j] = 1.0 / dt + std::abs(u[i][j]) * 3.0 / (2.0 * dx) + std::abs(v[i][j]) * 3.0 / (2.0 * dy)
                    + tempc * 2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)); //u3[i][j]の係数
                cx2[i][j] = -u[i][j] * 2.0 / (3.0 * dx) - std::abs(u[i][j]) / dx - tempc / (dx * dx); //u3[i-1][j]の係数
                cx3[i][j] = u[i][j] * 2.0 / (3.0 * dx) - std::abs(u[i][j]) / dx - tempc / (dx * dx); //u3[i+1][j]の係数
                cx4[i][j] = -v[i][j] * 2.0 / (3.0 * dy) - std::abs(v[i][j]) / dy - tempc / (dy * dy); //u3[i][j-1]の係数
                cx5[i][j] = v[i][j] * 2.0 / (3.0 * dy) - std::abs(v[i][j]) / dy - tempc / (dy * dy); //u3[i][j+1]の係数
                cx6[i][j] = u[i][j] / (12.0 * dx) + std::abs(u[i][j]) / (4.0 * dx); //u3[i-2][j]の係数
                cx7[i][j] = -u[i][j] / (12.0 * dx) + std::abs(u[i][j]) / (4.0 * dx); //u3[i+2][j]の係数
                cx8[i][j] = v[i][j] / (12.0 * dy) + std::abs(v[i][j]) / (4.0 * dy); //u3[i][j-2]の係数
                cx9[i][j] = -v[i][j] / (12.0 * dy) + std::abs(v[i][j]) / (4.0 * dy); //u3[i][j+2]の係数
            }
        }


        //BiCG法の初期設定
        #pragma omp for private(j) reduction(+:rr)
        for (i = 2; i < l + 1; i++) {
            for (j = 2; j < m + 1; j++) {
                double f = 0.0; //外力、もしくはソース項
                if (mode == 1) {
                    //vの場合のみ浮力を加える(y方向の浮力)
                    //f = beta * (tempera[i][j] - Tmin); //浮力と渦
                    f = beta * (tempera[i][j] - Tmin) + vortCon(i, j, false); //浮力と渦
                }
                else if (mode == 3) {
                    //燃料密度は時間経過に従って指数関数的に減少
                    f = -crho * rho[i][j];
                }
                else if (mode == 2) {
                    //燃料の密度に比例して熱を発生
                    f = rhot * rho[i][j];
                }
                else if (mode == 0) {
                    f = vortCon(i, j, true);
                }

                double xcon = ans[i][j] / dt + f; //cx1で割る前の右辺

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
        //y=Appを計算
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

        //修正係数αを計算
        alpha = rr / r0y;

        #pragma omp parallel
        {
            //sを計算
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    s[i][j] = resi[i][j] - alpha * y[i][j];
                }
            }

            //z=Asを計算
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

        //修正係数ωを計算
        omg = zs / zz;

        #pragma omp parallel
        {
            //次のステップの近似値を計算
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    ans3[i][j] += alpha * pp[i][j] + omg * s[i][j];
                }
            }

            //次のステップの近似に対する残差を計算
            #pragma omp for private(j) reduction(+:r2) reduction(+:nrr)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    resi[i][j] = s[i][j] - omg * z[i][j];
                    r2 += resi[i][j] * resi[i][j];
                    nrr += resia[i][j] * resi[i][j];
                }
            }
        }

        //残差が小さくなるかループ数が一定値以上になれば終了
        if (pow(r2, 0.5) < 1.0e-5 || cnt > l* m) {
            continueFlag = false;
        }

        //betaを計算
        beta = (alpha /omg) * (nrr / rr);

        //次のステップの方向ベクトルppを計算
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

    //境界条件を設定
    #pragma omp parallel
    {
        #pragma omp for
        for (i = 0; i < l + 3; i++) {
            ans3[i][0] = ans3[i][2]; //仮想セル
            ans3[i][1] = tempp; //下の面の温度Tmin
            ans3[i][m + 1] = tempp; //上の面の温度Tmin
            ans3[i][m + 2] = ans3[i][m]; //仮想セル
        }
        #pragma omp for
        for (j = 0; j < m + 3; j++) {
            ans3[0][j] = ans3[2][j]; //仮想セル
            ans3[1][j] = tempp; //横の面の温度Tmin
            ans3[l + 1][j] = tempp; //横の面の温度Tmin
            ans3[l + 2][j] = ans3[l][j]; //仮想セル
        }
        #pragma omp for
        for (j = 9 * m / 20; j < 11 * m / 20; j++) {
            ans3[1][j] = tempb; //吹き出し口の温度Tmax
        }
    }
}

//-----------------------------------------------------
//圧力項の計算(CG法)
//-----------------------------------------------------
void pressureByCG() {
    bool continueFlag = true;
    int cnt = 0;

    double h; //圧力のポアソン方程式の右辺

    double resi[l + 3][m + 3]; //点(i,j)における残差
    double pp[l + 3][m + 3]; //点(i,j)における修正方向ベクトル
    double y[l + 3][m + 3]; //App
    double alpha = 0.0; //解の修正係数α
    double beta = 0.0; //方向ベクトルpの修正係数β
    double rr = 0.0; //残差resiの内積
    double nrr = 0.0; //次のステップの残差resiの内積
    double ppy = 0.0; //ppとyの内積

    int i, j;

    //CG法の初期設定
    #pragma omp parallel
    {
        #pragma omp for private(j) reduction(+:rr)
        for (i = 2; i < l + 1; i++) {
            for (j = 2; j < m + 1; j++) {
                h = ((u3[i + 1][j] - u3[i - 1][j]) / dx + (v3[i][j + 1] - v3[i][j - 1]) / dy) / (2 * dt); //圧力のポアソン方程式の右辺
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

        //yを計算
        #pragma omp parallel for private(j) reduction(+:ppy)
        for ( i = 2; i < l + 1; i++) {
            for ( j = 2; j < m + 1; j++) {
                y[i][j] = cp1 * pp[i][j] + cp2 * pp[i - 1][j] + cp3 * pp[i + 1][j] + cp4 * pp[i][j - 1] + cp5 * pp[i][j + 1];
                ppy += pp[i][j] * y[i][j];
            }
        }

        //修正係数αを計算
        alpha = rr / ppy;

        #pragma omp parallel
        {
            //次のステップの近似値を計算
            #pragma omp for private(j)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    p[i][j] += alpha * pp[i][j];
                }
            }

            //次のステップの近似に対する残差を計算
            #pragma omp for private(j) reduction(+:nrr)
            for (i = 2; i < l + 1; i++) {
                for (j = 2; j < m + 1; j++) {
                    resi[i][j] -= alpha * y[i][j];
                    nrr += resi[i][j] * resi[i][j];
                }
            }
        }

        //残差が小さくなるかループ数が一定値以上になれば終了
        if (pow(nrr, 0.5) < 1.0e-5 || cnt > l* m) {
            continueFlag = false;
        }

        //betaを計算
        beta = nrr / rr;

        //次のステップの方向ベクトルppを計算
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

    //圧力勾配の計算
    #pragma omp parallel for private(j)
    for ( i = 1; i < l + 2; i++) {
        for ( j = 1; j < m + 2; j++) {
            px[i][j] = (p[i + 1][j] - p[i - 1][j]) / (2.0 * dx);
            py[i][j] = (p[i][j + 1] - p[i][j - 1]) / (2.0 * dy);
        }
    }
}

//-----------------------------------------------------
//速度の更新
//-----------------------------------------------------
void velo(bool isU) {
    double(*ans)[m + 3]; //u,vのどちらか
    double(*ans3)[m + 3]; //u3,v3のどちらか
    double(*ansp)[m + 3]; //px,pyのどちらか
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


        //速度の境界条件
        #pragma omp for
        for (int i = 0; i < l + 3; i++) {
            ans[i][0] = ans[i][2]; //仮想セル
            ans[i][1] = 0.0; //境界面の速度0
            ans[i][m + 1] = 0.0; //境界面の速度0
            ans[i][m + 2] = ans[i][m]; //仮想セル
        }
        #pragma omp for
        for (int j = 0; j < m + 3; j++) {
            ans[0][j] = ans[2][j]; //仮想セル
            ans[1][j] = 0.0; //境界面の速度0
            ans[l + 1][j] = 0.0; //境界面の速度0
            ans[l + 2][j] = ans[l][j]; //仮想セル
        }
        #pragma omp for
        for (int j = 9 * m / 20; j < 11 * m / 20; j++) {
            ans[1][j] = tempb; //吹き出し口のuはub、vはvb
        }
    }
}
