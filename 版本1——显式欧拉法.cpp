#include<stdio.h>
#include<math.h>
#include <graphics.h>
#include <conio.h>

#define G 6.67408     //��λ10e-11
#define dt 0.001    //������λʱ��


static int n=2;

double obj[2][8]={{0,0,0,0,0,0,591250,0},{384,0,0,102,0,0,42,0}};//�������ݵ�λ����ٶ����ݵ�λ����/s��������λ10e20
//(x����,y����,x�ٶ�,y�ٶ�,x���ٶ�,y���ٶ�,����,��־)
void Gravity(double obj[][8])
{
	int i,j;
	double RS;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i==j){continue;}
			RS=pow(obj[i][0]-obj[j][0],2)+pow(obj[i][1]-obj[j][1],2);
			obj[i][4]=obj[i][4]+G*obj[j][6]*(obj[j][0]-obj[i][0])/pow(RS,1.5);
			obj[i][5]=obj[i][5]+G*obj[j][6]*(obj[j][1]-obj[i][1])/pow(RS,1.5);
			}
	}
}

void CoVe(double obj[][8])
{
	int i=0;
	double t=0,x=0;
	for(i=0;i<n;i++){
		if(obj[i][7]==1){continue;}
		t=obj[i][2];
		obj[i][2]=obj[i][2]+obj[i][4]*dt;
		obj[i][0]=obj[i][0]+(obj[i][2]+t)*dt/2;
		x=obj[i][3];
		obj[i][3]=obj[i][3]+obj[i][5]*dt;
		obj[i][1]=obj[i][1]+(obj[i][3]+x)*dt/2;
		}
}

int main(){
	initgraph(1080, 720);
	while(1){
		Gravity(obj);
		CoVe(obj);
		setcolor(YELLOW);
		setfillcolor(GREEN);
		fillcircle(obj[0][0]/10+500, obj[0][1]/10+300, 5);
		setcolor(YELLOW);
		setfillcolor(BLUE);
		fillcircle(obj[1][0]/10+500, obj[1][1]/10+300, 2);
  // ��ʱ
	Sleep(5);   
}

return 0;
}