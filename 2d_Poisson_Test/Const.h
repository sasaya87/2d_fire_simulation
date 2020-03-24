#pragma once

const int Width = 640; //ウインドウ幅
const int Height = 480; //ウインドウ高さ
const float OrthWidth = 4.0f; //並行投影時の幅

const int l = 100; //x方向の分割数
const int m = 50; //y方向の分割数

const double dx = 0.02; //x方向のステップ幅
const double dy = 0.02; //y方向のステップ幅

const double dt = 0.001; //時間のステップ幅

const double Tmax = 5.0; //最高温度
const double Tmin = 1.0; //最低温度

const double pb = 1.0; //境界面の圧力

const double ub = 100.0; //uの吹き出し速度
const double vb = -40.0; //vの吹き出し速度

const double come = 300.0; //渦度強制法の係数

const double rhob = 1.0; //吹き出し時の燃料密度
const double crho = 20.0; //燃料の時間経過による減少の係数
const double arho = 0.1; //燃料の拡散にかかわる係数

const double betab = 300.0; //浮力の係数
const double a = 10.0; //熱の拡散にかかわる係数
const double rhot = 400.0; //燃料密度と発生熱量の比
const double re = 2500.0; //レイノルズ数