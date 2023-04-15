#include<stdio.h>
#include<math.h>
#include<float.h>

#include <graphics.h>
#include <conio.h>

#define G 6.67
#define p 10E-6
double t=0;

/**
struct Planet{
0	double vx;
1	double vy;
2	double x;
3	double y;
4	double mass;
5	int fix;
6	double K1;      X
7	double L1;
8	double K2;
9	double L2;
10	double K3;
11	double L3;
12	double K4;
13	double L4;
14	double K1;      Y
15	double L1;
16	double K2;
17	double L2;
18	double K3;
19	double L3;
20	double K4;
21	double L4;
};
struct Planet planet[4]={{123,32,43,2,5.965E24,0},{214,332,3844039,2,1023,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
**/

double planet[4][22]= {{-2,0,0,0,5708100,0},{172,-5,3,1400,59650,0},{100,0,0,1500,4500,},{98,0,0,1350,4841}};

double ffunc(double t,double x,double v){ //时间，原函数，一阶导
	return v;                
}

double gfunc(int i,int x_axis,int s, double h)
{                          //求二阶导，即加速度的函数；（对象在数组中的序号，计算方向x轴还是y轴，rk法公式中的预测步，步长）
	int j;
	double buff=0;
	double rs;
	for(j=0;j<4;j++)
	{
		if(i==j){continue;}
		if(x_axis==1)
		{//x轴向计算
			if(s==1){
				rs=pow(planet[i][2]-planet[j][2],2)+pow(planet[i][3]-planet[j][3],2);
				if(rs<p){break;}
                buff+=G*planet[j][4]*(planet[j][2]-planet[i][2])/pow(rs,1.5);
			}
			if(s==2){
				rs=pow(planet[i][2]+h*planet[i][6]/2-planet[j][2]-h*planet[j][6]/2,2)+pow(planet[i][3]+h*planet[i][14]/2-planet[j][3]-h*planet[j][14]/2,2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][2]+h*planet[j][6]/2-planet[i][2]-h*planet[i][6]/2)/pow(rs,1.5);		
			}
			if(s==3){
				rs=pow(planet[i][2]+h*planet[i][8]/2-planet[j][2]-h*planet[j][8]/2,2)+pow(planet[i][3]+h*planet[i][16]/2-planet[j][3]-h*planet[j][16]/2,2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][2]+h*planet[j][8]/2-planet[i][2]-h*planet[i][8]/2)/pow(rs,1.5);
			}
			if(s==4){
				rs=pow(planet[i][2]+h*planet[i][10]-planet[j][2]-h*planet[j][10],2)+pow(planet[i][3]+h*planet[i][18]-planet[j][3]-h*planet[j][18],2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][2]+h*planet[j][10]-planet[i][2]-h*planet[i][10])/pow(rs,1.5);
			}
		
		}
		if(x_axis==0)
		{//y轴向计算
			if(s==1){
				rs=pow(planet[i][3]-planet[j][3],2)+pow(planet[i][2]-planet[j][2],2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][3]-planet[i][3],2)/pow(rs,1.5);
			}
			if(s==2){
				rs=pow(planet[i][3]+h*planet[i][14]/2-planet[j][3]-h*planet[j][14]/2,2)+pow(planet[i][2]+h*planet[i][6]/2-planet[j][2]-h*planet[j][6]/2,2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][3]+h*planet[j][14]/2-planet[i][3]-h*planet[i][14]/2)/pow(rs,1.5);	
			}
			if(s==3){
				rs=pow(planet[i][2]+h*planet[i][8]/2-planet[j][2]-h*planet[j][8]/2,2)+pow(planet[i][3]+h*planet[i][16]/2-planet[j][3]-h*planet[j][16]/2,2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][3]+h*planet[j][16]/2-planet[i][3]-h*planet[i][16]/2)/pow(rs,1.5);
			}
			if(s==4){
				rs=pow(planet[i][2]+h*planet[i][10]-planet[j][2]-h*planet[j][10],2)+pow(planet[i][3]+h*planet[i][18]-planet[j][3]-h*planet[j][18],2);
				if(rs<p){break;}
				buff+=G*planet[j][4]*(planet[j][3]+h*planet[j][18]-planet[i][3]-h*planet[i][18])/pow(rs,1.5);
			}
		}
	}
	return buff;
}

double rk4(double h){
	int i=0;
	for(i=0;i<4;i++){
		planet[i][6]=ffunc(t,planet[i][2],planet[i][0]);
		planet[i][7]=gfunc(i,1,1,h);
		planet[i][14]=ffunc(t,planet[i][3],planet[i][1]);
		planet[i][15]=gfunc(i,0,1,h);}
	for(i=0;i<4;i++){
		planet[i][8]=ffunc(t+h/2,planet[i][2]+h*planet[i][6],planet[i][0]+h*planet[i][7]/2);
		planet[i][9]=gfunc(i,1,2,h);
		planet[i][16]=ffunc(t+h/2,planet[i][3]+h*planet[i][14],planet[i][1]+h*planet[i][15]/2);
		planet[i][17]=gfunc(i,0,2,h);}
	for(i=0;i<4;i++){
		planet[i][10]=ffunc(t+h/2,planet[i][2]+h*planet[i][8]/2,planet[i][0]+h*planet[i][9]/2);
		planet[i][11]=gfunc(i,1,3,h);
		planet[i][18]=ffunc(t+h/2,planet[i][3]+h*planet[i][16]/2,planet[i][1]+h*planet[i][17]/2);
		planet[i][19]=gfunc(i,0,3,h);}
	for(i=0;i<4;i++){
		planet[i][12]=ffunc(t+h,planet[i][2]+h*planet[i][10],planet[i][0]+h*planet[i][11]);
		planet[i][13]=gfunc(i,1,4,h);
		planet[i][20]=ffunc(t+h,planet[i][3]+h*planet[i][18],planet[i][1]+h*planet[i][19]);
		planet[i][21]=gfunc(i,0,4,h);}
	for(i=0;i<4;i++){
		planet[i][2]+=(planet[i][6]+2*planet[i][8]+2*planet[i][10]+planet[i][12])*h/6;	
		planet[i][0]+=(planet[i][7]+2*planet[i][9]+2*planet[i][11]+planet[i][13])*h/6;
		planet[i][3]+=(planet[i][14]+2*planet[i][16]+2*planet[i][18]+planet[i][20])*h/6;
		planet[i][1]+=(planet[i][15]+2*planet[i][17]+2*planet[i][19]+planet[i][21])*h/6;
	}
	return 0;
}

void drawpath(int i,char color){
	
}

int main(){

		initgraph(1080,720);	
	while(1){
		rk4(0.01);
		setcolor(GREEN);
		setfillcolor(GREEN);
		fillcircle((int)planet[0][2]/10+500, (int)planet[0][3]/10+300, 1);
		setcolor(BLUE);
		setfillcolor(BLUE);
		fillcircle((int)planet[1][2]/10+500, (int)planet[1][3]/10+300, 2);
		setcolor(YELLOW);
		setfillcolor(YELLOW);
		fillcircle((int)planet[2][2]/10+500, (int)planet[2][3]/10+300, 1);
		setcolor(RED);
		setfillcolor(RED);
		fillcircle((int)planet[3][2]/10+500, (int)planet[3][3]/10+300, 1);
		Sleep(1);
	}
	return 0;
}