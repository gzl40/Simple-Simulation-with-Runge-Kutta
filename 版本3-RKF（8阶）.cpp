#include<stdio.h>
#include<math.h>
#include<float.h>

#include <graphics.h>
#include <conio.h>
#define n 4
#define G 6.67
#define p1 10E-4//浮点数等于零的判定精度
#define p2 10E-6//自适应步长的判定值


/**
struct Planet{
0	double vx;
1	double vy;
2	double x;
3	double y;
4	double mass;
5	int fix;

K1     0//a
K2
K3
K4
K5
K6
K7
K8
K9
K10
K11
K12
K13
L1      13//v
L2
L3
L4
L5
L6
L7
L8     20
L9     21
L10    22
L11    23
L12    24
L13    25

};
struct Planet planet[4]={{123,32,43,2,5.965E24,0},{214,332,3844039,2,1023,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
**/

double planet[n][6]={{0,0,0,0,5708100,0},{172,-5,3,1400,59650,0},{100,0,0,1500,4500,0},{98,0,0,1350,4841,0}};
double xc[n][26]={{},{},{},{}};
double yc[n][26]={{},{},{},{}};
double er[2*n]={};


double rk8(double h){
	int i=0;      //rk步数
	int j=0;        //万有引力定律另一个对象
	double x;   //x向距离
	double y;
	double buff;//储存合加速度
	double rs;    //距离的平方
	double ai;    //系数
	for(i=0;i<n;i++){//1
		ai=0;
		xc[i][13]=planet[i][0];
		yc[i][13]=planet[i][1];
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]-planet[i][2];
			y=planet[j][3]-planet[i][3];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5); 
			xc[i][0]+=x*buff;
			yc[i][0]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//2
		ai=2/27.0;
		xc[i][14]=planet[i][0]+ai*xc[i][0]*h;
		yc[i][14]=planet[i][1]+ai*yc[i][0]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(2/27.0)*h*xc[j][13]-planet[i][2]-(2/27.0)*h*xc[i][13];
			y=planet[j][3]+(2/27.0)*h*yc[j][13]-planet[i][3]-(2/27.0)*h*yc[i][13];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);
			xc[i][1]+=x*buff;
			yc[i][1]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//3
   		ai=1/9.0;
		xc[i][15]=planet[i][0]+ai*xc[i][1]*h;
		yc[i][15]=planet[i][1]+ai*yc[i][1]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(1/36.0)*h*xc[j][13]+(1/12.0)*h*xc[j][14]-planet[i][2]-(1/36.0)*h*xc[i][13]-(1/12.0)*h*xc[i][14];
			y=planet[j][3]+(1/36.0)*h*yc[j][13]+(1/12.0)*h*yc[j][14]-planet[i][3]-(1/36.0)*h*yc[i][13]-(1/12.0)*h*yc[i][14];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][2]+=x*buff;
			yc[i][2]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//4
		ai=1/6.0;
		xc[i][16]=planet[i][0]+ai*xc[i][2]*h;
		yc[i][16]=planet[i][1]+ai*yc[i][2]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(1/24.0)*h*xc[j][13]+(1/8.0)*h*xc[j][15]-planet[i][2]-(1/24.0)*h*xc[i][13]-(1/8.0)*h*xc[i][15];
			y=planet[j][3]+(1/24.0)*h*yc[j][13]+(1/8.0)*h*yc[j][15]-planet[i][3]-(1/24.0)*h*yc[i][13]-(1/8.0)*h*yc[i][15];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][3]+=x*buff;
			yc[i][3]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//5
		ai=5/12.0;
		xc[i][17]=planet[i][0]+ai*xc[i][3]*h;
		yc[i][17]=planet[i][1]+ai*yc[i][3]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(5/12.0)*h*xc[j][13]+(-25/16.0)*h*xc[j][15]+(25/16.0)*h*xc[j][16]-planet[i][2]-(5/12.0)*h*xc[i][13]-(-25/16.0)*h*xc[i][15]-(25/16.0)*h*xc[i][16];
			y=planet[j][3]+(5/12.0)*h*yc[j][13]+(-25/16.0)*h*yc[j][15]+(25/16.0)*h*yc[j][16]-planet[i][3]-(5/12.0)*h*yc[i][13]-(-25/16.0)*h*yc[i][15]-(25/16.0)*h*yc[i][16];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][4]+=x*buff;
			yc[i][4]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//6
		ai=1/2.0;
		xc[i][18]=planet[i][0]+ai*xc[i][4]*h;
		yc[i][18]=planet[i][1]+ai*yc[i][4]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(1/20.0)*h*xc[j][13]+(1/4.0)*h*xc[j][16]+(1/5.0)*h*xc[j][17]-planet[i][2]-(1/20.0)*h*xc[i][13]-(1/4.0)*h*xc[j][16]-(1/5.0)*h*xc[i][17];
			y=planet[j][3]+(1/20.0)*h*yc[j][13]+(1/4.0)*h*yc[j][16]+(1/5.0)*h*yc[j][17]-planet[i][3]-(1/20.0)*h*yc[i][13]-(1/4.0)*h*yc[j][16]-(1/5.0)*h*yc[i][17];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][5]+=x*buff;
			yc[i][5]+=y*buff;            
		}
	}
	for(i=0;i<n;i++){//7
		ai=5/6.0;
		xc[i][19]=planet[i][0]+ai*xc[i][5]*h;
		yc[i][19]=planet[i][1]+ai*yc[i][5]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(-25/108.0)*h*xc[j][13]+(125/108.0)*h*xc[j][16]+(-65/27.0)*h*xc[j][17]+(125/54.0)*h*xc[j][18]-planet[i][2]-(-25/108.0)*h*xc[i][13]-(125/108.0)*h*xc[j][16]-(-65/27.0)*h*xc[i][17]-(125/54.0)*h*xc[i][18];
			y=planet[j][3]+(-25/108.0)*h*yc[j][13]+(125/108.0)*h*yc[j][16]+(-65/27.0)*h*yc[j][17]+(125/54.0)*h*yc[j][18]-planet[i][3]-(-25/108.0)*h*yc[i][13]-(125/108.0)*h*yc[j][16]-(-65/27.0)*h*yc[i][17]-(125/54.0)*h*yc[i][18];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][6]+=x*buff;
			yc[i][6]+=y*buff;            
		}
	}
	for(i=0;i<n;i++)//8
	{
        ai=1/6.0;
	    xc[i][20]=planet[i][0]+ai*xc[i][6]*h;
	    yc[i][20]=planet[i][1]+ai*yc[i][6]*h;
	    for(j=0;j<n;j++)
		{
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(31/300.0)*h*xc[j][13]+(61/225.0)*h*xc[j][17]+(-2/9.0)*h*xc[j][18]+(13/900.0)*h*xc[j][19]-
				(planet[i][2]+(31/300.0)*h*xc[i][13]+(61/225.0)*h*xc[i][17]+(-2/9.0)*h*xc[i][18]+(13/900.0)*h*xc[i][19]);
			y=planet[j][3]+(31/300.0)*h*yc[j][13]+(61/225.0)*h*yc[j][17]+(-2/9.0)*h*yc[j][18]+(13/900.0)*h*yc[j][19]-
				(planet[i][3]+(31/300.0)*h*yc[i][13]+(61/225.0)*h*yc[i][17]+(-2/9.0)*h*yc[i][18]+(13/900.0)*h*yc[i][19]);
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][7]+=x*buff;
			yc[i][7]+=y*buff;         		
		}
	}
	for(i=0;i<n;i++)//9
	{
		ai=2/3.0;
	    xc[i][21]=planet[i][0]+ai*xc[i][7]*h;
	    yc[i][21]=planet[i][1]+ai*yc[i][7]*h;
	    for(j=0;j<n;j++)
		{
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(2)*h*xc[j][13]+(-53/6.0)*h*xc[j][16]+(704/45.0)*h*xc[j][17]+(-107/9.0)*h*xc[j][18]+(67/90.0)*h*xc[j][19]+(3)*h*xc[j][20]-
				(planet[i][2]+(2)*h*xc[i][13]+(-53/6.0)*h*xc[i][16]+(704/45.0)*h*xc[i][17]+(-107/9.0)*h*xc[i][18]+(67/90.0)*h*xc[i][19]+(3)*h*xc[i][20]);
			y=planet[j][3]+(2)*h*yc[j][13]+(-53/6.0)*h*yc[j][16]+(704/45.0)*h*yc[j][17]+(-107/9.0)*h*yc[j][18]+(67/90.0)*h*yc[j][19]+(3)*h*yc[j][20]-
				(planet[i][3]+(2)*h*yc[i][13]+(-53/6.0)*h*yc[i][16]+(704/45.0)*h*yc[i][17]+(-107/9.0)*h*yc[i][18]+(67/90.0)*h*yc[i][19]+(3)*h*yc[i][20]);
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][8]+=x*buff;
			yc[i][8]+=y*buff;         		
		}
	}
	for(i=0;i<n;i++)//10
	{
		ai=1/3.0;
	    xc[i][22]=planet[i][0]+ai*xc[i][8]*h;
	    yc[i][22]=planet[i][1]+ai*yc[i][8]*h;
	    for(j=0;j<n;j++)
		{
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(-91/108.0)*h*xc[j][13]+(23/108.0)*h*xc[j][16]+(-976/135.0)*h*xc[j][17]+(311/54.0)*h*xc[j][18]+(-19/60.0)*h*xc[j][19]+(17/6.0)*h*xc[j][20]+(-1/12.0)*h*xc[j][21]-
				(planet[i][2]+(-91/108.0)*h*xc[i][13]+(23/108.0)*h*xc[i][16]+(-976/135.0)*h*xc[i][17]+(311/54.0)*h*xc[i][18]+(-19/60.0)*h*xc[i][19]+(17/6.0)*h*xc[i][20]+(-1/12.0)*h*xc[i][21]);
			y=planet[j][3]+(-91/108.0)*h*yc[j][13]+(23/108.0)*h*yc[j][16]+(-976/135.0)*h*yc[j][17]+(311/54.0)*h*yc[j][18]+(-19/60.0)*h*yc[j][19]+(17/6.0)*h*yc[j][20]+(-1/12.0)*h*yc[j][21]-
				(planet[i][3]+(-91/108.0)*h*yc[i][13]+(23/108.0)*h*yc[i][16]+(-976/135.0)*h*yc[i][17]+(311/54.0)*h*yc[i][18]+(-19/60.0)*h*yc[i][19]+(17/6.0)*h*yc[i][20]+(-1/12.0)*h*yc[i][21]);
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][9]+=x*buff;
			yc[i][9]+=y*buff;         		
		}
	}
	for(i=0;i<n;i++)//11
    {
		ai=1;
	    xc[i][23]=planet[i][0]+ai*xc[i][9]*h;
	    yc[i][23]=planet[i][1]+ai*yc[i][9]*h;
	    for(j=0;j<n;j++)
		{
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(2383/4100.0)*h*xc[j][13]+(-341/164.0)*h*xc[j][16]+(4496/1025.0)*h*xc[j][17]+(-301/82.0)*h*xc[j][18]+(2133/4100.0)*h*xc[j][19]+(45.0/82.0)*h*xc[j][20]+(45.0/164.0)*h*xc[j][21]+(18/41.0)*h*xc[j][22]-
				(planet[i][2]+(2383/4100.0)*h*xc[i][13]+(-341/164.0)*h*xc[i][16]+(4496/1025.0)*h*xc[i][17]+(-301/82.0)*h*xc[i][18]+(2133/4100.0)*h*xc[i][19]+(45.0/82.0)*h*xc[i][20]+(45.0/164.0)*h*xc[i][21]+(18/41.0)*h*xc[i][22]);
			y=planet[j][3]+(2383/4100.0)*h*yc[j][13]+(-341/164.0)*h*yc[j][16]+(4496/1025.0)*h*yc[j][17]+(-301/82.0)*h*yc[j][18]+(2133/4100.0)*h*yc[j][19]+(45.0/82.0)*h*yc[j][20]+(45.0/164.0)*h*yc[j][21]+(18/41.0)*h*yc[j][22]-
				(planet[i][3]+(2383/4100.0)*h*yc[i][13]+(-341/164.0)*h*yc[i][16]+(4496/1025.0)*h*yc[i][17]+(-301/82.0)*h*yc[i][18]+(2133/4100.0)*h*yc[i][19]+(45.0/82.0)*h*yc[i][20]+(45.0/164.0)*h*yc[i][21]+(18/41.0)*h*yc[i][22]);
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][10]+=x*buff;
			yc[i][10]+=y*buff;         		
		}
	}
	for(i=0;i<n;i++){//12
		ai=0;
		xc[i][24]=planet[i][0]+ai*xc[i][10]*h;
		yc[i][24]=planet[i][1]+ai*yc[i][10]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(3/205.0)*h*xc[j][13]+(-6/41.0)*h*xc[j][18]+(-3/205.0)*h*xc[j][19]+(-3/41.0)*h*xc[j][20]+(3/41.0)*h*xc[j][21]+(6/41.0)*h*xc[j][22]-
				planet[i][2]-(3/205.0)*h*xc[i][13]-(-6/41.0)*h*xc[i][18]-(-3/205.0)*h*xc[i][19]-(-3/41.0)*h*xc[i][20]-(3/41.0)*h*xc[i][21]-(6/41.0)*h*xc[i][22];
			y=planet[j][3]+(3/205.0)*h*yc[j][13]+(-6/41.0)*h*yc[j][18]+(-3/205.0)*h*yc[j][19]+(-3/41.0)*h*yc[j][20]+(3/41.0)*h*yc[j][21]+(6/41.0)*h*yc[j][22]-
				planet[i][3]-(3/205.0)*h*yc[i][13]-(-6/41.0)*h*yc[i][18]-(-3/205.0)*h*yc[i][19]-(-3/41.0)*h*yc[i][20]-(3/41.0)*h*yc[i][21]-(6/41.0)*h*yc[i][22];
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][11]+=x*buff;
			yc[i][11]+=y*buff;         		
		}
	}
	for(i=0;i<n;i++){//13
		ai=1;
		xc[i][25]=planet[i][0]+ai*xc[i][11]*h;
		yc[i][25]=planet[i][1]+ai*yc[i][11]*h;
		for(j=0;j<n;j++){
			if(i==j){continue;}//跳过自己
			x=planet[j][2]+(-1777/4100.0)*h*xc[j][13]+(341.0/164.0)*h*xc[j][16]+(4496/1025.0)*h*xc[j][17]+(-289/82.0)*h*xc[j][18]+(2193/4100.0)*h*xc[j][19]+(51/82.0)*h*xc[j][20]+(33/164.0)*h*xc[j][21]+(12.0/41.0)*h*xc[j][22]+xc[j][24]*h-
				planet[i][2]-(-1777/4100.0)*h*xc[i][13]-(341.0/164.0)*h*xc[i][16]-(4496/1025.0)*h*xc[i][17]-(-289/82.0)*h*xc[i][18]-(2193/4100.0)*h*xc[i][19]-(51/82.0)*h*xc[i][20]-(33/164.0)*h*xc[i][21]-(12.0/41.0)*h*xc[i][22]-xc[i][24]*h;
			y=planet[j][3]+(-1777/4100.0)*h*yc[j][13]+(341.0/164.0)*h*yc[j][16]+(4496/1025.0)*h*yc[j][17]+(-289/82.0)*h*yc[j][18]+(2193/4100.0)*h*yc[j][19]+(51/82.0)*h*yc[j][20]+(33/164.0)*h*yc[j][21]+(12.0/41.0)*h*yc[j][22]+yc[j][24]*h-
				planet[i][3]-(-1777/4100.0)*h*yc[i][13]-(341.0/164.0)*h*yc[i][16]-(4496/1025.0)*h*yc[i][17]-(-289/82.0)*h*yc[i][18]-(2193/4100.0)*h*yc[i][19]-(51/82.0)*h*yc[i][20]-(33/164.0)*h*yc[i][21]-(12.0/41.0)*h*yc[i][22]-yc[i][24]*h;
			rs=pow(x,2)+pow(y,2);
			if(rs<p1){break;}
			buff=G*planet[j][4]/pow(rs,1.5);        
			xc[i][12]+=x*buff;
			yc[i][12]+=y*buff;      
		}
	}
	for(i=0;i<n;i++){
		planet[i][2]+=((0)*h*xc[i][13]+0*xc[i][14]+0*xc[i][15]+0*xc[i][16]+0*xc[i][17]+(34/105.0)*h*xc[i][18]+(9/35.0)*h*xc[i][19]+(9/35.0)*h*xc[i][20]+(9/280.0)*h*xc[i][21]+(9/280.0)*h*xc[i][22]+(0)*h*xc[i][23]+(41.0/840.0)*h*xc[i][24]+(41.0/840.0)*h*xc[i][25]);
		planet[i][3]+=((0)*h*yc[i][13]+0*yc[i][14]+0*yc[i][15]+0*yc[i][16]+0*yc[i][17]+(34/105.0)*h*yc[i][18]+(9/35.0)*h*yc[i][19]+(9/35.0)*h*yc[i][20]+(9/280.0)*h*yc[i][21]+(9/280.0)*h*yc[i][22]+(0)*h*yc[i][23]+(41.0/840.0)*h*yc[i][24]+(41.0/840.0)*h*yc[i][25]);
		planet[i][0]+=((0)*h*xc[i][0]+0*xc[i][1]+0*xc[i][2]+0*xc[i][3]+0*xc[i][4]+(34/105.0)*h*xc[i][5]+(9/35.0)*h*xc[i][6]+(9/35.0)*h*xc[i][7]+(9/280.0)*h*xc[i][8]+(9/280.0)*h*xc[i][9]+(0)*h*xc[i][10]+(41.0/840.0)*h*xc[i][11]+(41.0/840.0)*h*xc[i][12]);
		planet[i][1]+=((0)*h*yc[i][0]+0*yc[i][1]+0*yc[i][2]+0*yc[i][3]+0*yc[i][4]+(34/105.0)*h*yc[i][5]+(9/35.0)*h*yc[i][6]+(9/35.0)*h*yc[i][7]+(9/280.0)*h*yc[i][8]+(9/280.0)*h*yc[i][9]+(0)*h*yc[i][10]+(41.0/840.0)*h*yc[i][11]+(41.0/840.0)*h*yc[i][12]);
		er[2*i]=xc[i][0]+xc[i][10]-xc[i][11]-xc[i][12];
		er[2*i+1]=yc[i][0]+yc[i][10]-yc[i][11]-yc[i][12];
		xc[i][0]=xc[i][1]=xc[i][2]=xc[i][3]=xc[i][4]=xc[i][5]=xc[i][6]=xc[i][7]=xc[i][8]=xc[i][9]=xc[i][10]=xc[i][11]=xc[i][12]=0;
		yc[i][0]=yc[i][1]=yc[i][2]=yc[i][3]=yc[i][4]=yc[i][5]=yc[i][6]=yc[i][7]=yc[i][8]=yc[i][9]=yc[i][10]=yc[i][11]=yc[i][12]=0;
		xc[i][13]=xc[i][14]=xc[i][15]=xc[i][16]=xc[i][17]=xc[i][18]=xc[i][19]=xc[i][20]=xc[i][21]=xc[i][22]=xc[i][23]=xc[i][24]=xc[i][25]=0;
		yc[i][13]=yc[i][14]=yc[i][15]=yc[i][16]=yc[i][17]=yc[i][18]=yc[i][19]=yc[i][20]=yc[i][21]=yc[i][22]=yc[i][23]=yc[i][24]=yc[i][25]=0;
	}
	return 0;
}

double erC(){
	double er[2*n]={},avg=0;
	int i;
	for(i=0;i<2*n;i++){avg+=er[i];}
	avg=avg/(2*n);
	return avg;
}

double AutoFoto(double h){
	double fo=h,k=1;
	while(1)
	{
		rk8(k*fo);
		if(erC()<0.1)
		{
			if(0.95-k<p2){break;}
			else{k=(0.5*(0.95/k-1)+1)*k;}
		}
		else
		{
			if(k-0.1<p2){break;}
			else{k=(0.5*((k/0.1)-1)+1)*0.1;}
		}
	}
	return 0;
}

void button(int x,int y,int w,int h,TCHAR * text)
 {
    setbkmode(TRANSPARENT);
    setfillcolor(GREEN);
    fillroundrect(x,y,x+w,y+h,10,10);
    // 输出字符串（MBCS 字符集）
    TCHAR s1[] = L"黑体";
    settextstyle(15,0,s1);
    TCHAR s[50] = L"hello";
    
    int tx = x + (w - textwidth(text)) / 2;
    int ty = y + (h - textheight(text)) / 2;

    outtextxy(tx, ty, text);
}

void drawpath(int i,char color){
}
void AllB(){
		TCHAR s[50] = L"Pause";
		TCHAR s1[50] = L"Pace down";
		TCHAR s2[50] = L"Pace up";
		TCHAR s3[50] = L"Size down";
		TCHAR s4[50] = L"Size up";
		TCHAR s5[50] = L"Shift L";
		TCHAR s6[50] = L"Shift U";
		TCHAR s7[50] = L"Shift D";
		TCHAR s8[50] = L"Shift R";
		button(10,10,85,25,s1);
		button(95, 10, 85, 25, s);
		button(180,10,85,25,s2);
		button(10, 39,127.5, 25, s3);
		button(137.5,39,127.5, 25, s4);
		button(10,70,85, 50, s5);
		button(95,70,85, 25, s6);
		button(95,95,85, 25, s7);
		button(180,70,85, 50, s8);
		settextcolor(RED);//设置文字的颜色，这里为红色

		//outtextxy(50,50,L"123");//前两个数字为输入文字的坐标，双引号中填输出的东西
}
int main(){
	int tag=1,p=10,rate=10,i;
	int xoff=500;
	int yoff=300;
		initgraph(1080,720);	
		cleardevice();
		ExMessage msg;
		AllB();
	while(1){
     if (peekmessage(&msg, EM_MOUSE)) {
         switch (msg.message)
         {
         case WM_LBUTTONDOWN:
             if (msg.x >= 10 && msg.x <= 10+85 && msg.y >= 10 && msg.y <= 10+25)
             {
				 if(p<2){break;}
				p-=2;
             }
             if (msg.x >= 95 && msg.x <= 95+85 && msg.y >= 10 && msg.y <= 10+25)
             {
				if(tag==1){tag=0;}
				else tag=1;	
             }    
			 if (msg.x >= 180 && msg.x <= 180+85 && msg.y >= 10 && msg.y <= 10+25)
             {
				p+=10;
             }
			 if (msg.x >= 10 && msg.x <= 10+85 && msg.y >= 37 && msg.y <= 37+25)
             {
				rate+=5;
				initgraph(1080,720);
				AllB();
             }
			 if (msg.x >= 137.5 && msg.x <= 137.5+85 && msg.y >= 37 && msg.y <= 37+25)
             {
				 if(rate<6){break;}
				 rate-=5;
				initgraph(1080,720);
				AllB();
             }
			 if (msg.x >= 10 && msg.x <= 10+85 && msg.y >= 70 && msg.y <= 70+50)
             {
				xoff+=25;
				initgraph(1080,720);
				AllB();
             }
			 if (msg.x >= 95 && msg.x <= 95+85 && msg.y >= 70 && msg.y <= 70+25)
             {
				yoff+=25;
				initgraph(1080,720);
				AllB();
             }
			 if (msg.x >= 95 && msg.x <= 95+85 && msg.y >= 95 && msg.y <= 95+25)
             {
				yoff-=25;
				initgraph(1080,720);
				AllB();
             }
			 if (msg.x >= 180 && msg.x <= 180+85 && msg.y >= 70 && msg.y <= 70+50)
             {
				xoff-=25;
				initgraph(1080,720);
				AllB();
             }
             break;
         default:
             break;
         }
    }
	 if(tag==1){for(i=0;i<p;i++){AutoFoto(0.0001);}}
		setcolor(GREEN);
		setfillcolor(GREEN);
		fillcircle((int)planet[0][2]/rate+xoff, (int)planet[0][3]/rate+yoff, 1);
		setcolor(WHITE);
		setfillcolor(WHITE);
		fillcircle((int)planet[1][2]/rate+xoff, (int)planet[1][3]/rate+yoff, 1);
		setcolor(YELLOW);
		setfillcolor(YELLOW);
		fillcircle((int)planet[2][2]/rate+xoff, (int)planet[2][3]/rate+yoff, 1);
		setcolor(RED);
		setfillcolor(RED);
		fillcircle((int)planet[3][2]/rate+xoff, (int)planet[3][3]/rate+yoff, 1);
		}
	return 0;
}