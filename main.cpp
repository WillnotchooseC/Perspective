#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cstdlib> 
#include <cstdio> 
#include <cmath> 
#include <fstream> 
#include <vector> 
#include <iostream> 
#include <cassert> 
using namespace std;
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) < (y) ? (y) : (x))
//#define M_PI 3.141592653589793 
#define INTERVAL 0.2f
#define IndexVect(i,j,k) ((k)*128*128 + (j)*128 + (i))
//#include <sys/times.h>
//using namespace std;
//#define IX(size, i, j, k) ( i * j * k * size + i * j * size + i ) 
static float d = 2.0f;
int xdim, ydim, zdim;
typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} RGBDATA;
static const float XR = 0.0f, YR = 0.0f, ZR = 0.0f;



void WritePPM(int, int, RGBDATA *, FILE*);
void ReadRawiv(int* , int*, int*, float**, char *);

/*template<typename T>
class Vec3
{
public:
	T x, y, z;
	Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
	Vec3(T xx) : x(xx), y(xx), z(xx) {}
	Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
	Vec3& normalize()
	{
		T nor2 = length2();
		if (nor2 > 0) {
			T invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}
	Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
	Vec3<T> operator / (const T &f) const { return Vec3<T>(x / f, y / f, z / f); }
	Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
	T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
	Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
	Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
	Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
	Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
	T length2() const { return x * x + y * y + z * z; }
	T length() const { return sqrt(length2()); }
	friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
	{
		os << "[" << v.x << " " << v.y << " " << v.z << "]";
		return os;
	}
};*/

//typedef Vec3<float> Vec3f;
struct Vec3f
{
	float x, y, z;
	Vec3f operator+(const Vec3f &v)
	{
		return{ x + v.x, y + v.y, z + v.z };
	}
	Vec3f operator-(const Vec3f &v)
	{
		return{ x - v.x, y - v.y, z - v.z };
	}
	Vec3f operator*(const float val)
	{
		return{ x*val, y*val, z*val };
	}

	Vec3f& operator *= (const Vec3f &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vec3f operator/(const float val)
	{
		return{ x / val, y / val, z / val };
	}
	float dot(const Vec3f &v)
	{
		return x*v.x + y*v.y + z*v.z;
	}
	float length2() const { return x * x + y * y + z * z; }
	Vec3f& normalize()
	{
		float nor2 = length2();
		if (nor2 > 0) {
			float invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}
};
struct matrix
{
	Vec3f r1, r2, r3;
	matrix operator*(const matrix &m)
	{
		return{
			{ r1.dot({m.r1.x, m.r2.x, m.r3.x}),
			r1.dot({m.r1.y, m.r2.y, m.r3.y }),
			r1.dot({m.r1.z, m.r2.z, m.r3.z })
			},
			{ r2.dot({m.r1.x, m.r2.x, m.r3.x }),
			r2.dot({m.r1.y, m.r2.y, m.r3.y }),
			r2.dot({m.r1.z, m.r2.z, m.r3.z })
			},
			{ r3.dot({m.r1.x, m.r2.x, m.r3.x }),
			r3.dot({m.r1.y, m.r2.y, m.r3.y }),
			r3.dot({m.r1.z, m.r2.z, m.r3.z })
			},
		};
	}
};

Vec3f mul(Vec3f &vector, matrix &m) {
	return{ vector.dot({m.r1.x , m.r2.x , m.r3.x}),
			vector.dot({m.r1.y , m.r2.y , m.r3.y}),
			vector.dot({m.r1.z , m.r2.z , m.r3.z}) };
}


/*Vec3f trace(const Vec3f &rayorig, const Vec3f &raydir, float *volume, float *alpha, const Vec3f *color, float& Irange) {
	Vec3f lookDirection = raydir -rayorig;
	lookDirection.normalize();
	//printf("in the trace ");
	std::cout << "in the trace";
	float newalpha = 0.0f;
	float currentValue = 0.0f;
	int time = 1;
	Vec3f currentColor = Vec3f(0);
	Vec3f finalColor = Vec3f(0);
	const float stepD = 10.0f;
	std::cout<<volume[130];
	do {
		
		lookDirection = lookDirection*stepD*time;
		int x0 = (int)(lookDirection.x / 1.768193);
		int x1 = x0 + 1;
		int y0 = (int)(lookDirection.y / 1.768193);
		int y1 = y0 + 1;
		int z0 = (int)(lookDirection.z / 1.444882);
		int z1 = z0 - 1;
		

		currentValue = volume[IX(127,x0, y0, z0)] * (1 - (lookDirection.x / 224.560547))*(1 - (lookDirection.y / 224.560547))
			*(1 - (lookDirection.z / 183.500000)) + volume[IX(127, x1, y0, z0)] * (lookDirection.x / 224.560547)*
			(1 - (lookDirection.y / 224.560547))*(1 - (lookDirection.z / 183.500000)) + volume[IX(127, x0, y1, z0)] *
			(1 - (lookDirection.x / 224.560547))*(lookDirection.y / 224.560547)*(1 - (lookDirection.z / 183.500000)) +
			volume[IX(127, x0, y0, z1)] * (1 - (lookDirection.x / 224.560547))*(1 - (lookDirection.y / 224.560547))
			*(lookDirection.z / 183.500000) + volume[IX(127, x1, y0, z1)] * (lookDirection.x / 224.560547)
			*(1 - (lookDirection.y / 224.560547))*(lookDirection.z / 183.500000) + volume[IX(127, x0, y1, z1)] *
			(1 - (lookDirection.x / 224.560547))*(lookDirection.y / 224.560547)*(lookDirection.z / 183.500000) +
			volume[IX(127, x1, y1, z0)] * (lookDirection.x / 224.560547)*(lookDirection.y / 224.560547)*(1 - (lookDirection.z / 183.500000)) +
			volume[IX(127, x1, y1, z1)] * (lookDirection.x / 224.560547)*(lookDirection.y / 224.560547)*(lookDirection.z / 183.500000);
		//printf("this is currentValue \n");
		//printf("%d %d %d %d %d %d\n", x0,x1,y0,y1,z0,z1);
		currentColor = color[(int)((currentValue / (Irange+0.2)) * 5)]+ currentColor*(1 - alpha[(int)((currentValue / (Irange + 0.2)) * 5)]);
		newalpha = alpha[(int)((currentValue / (Irange + 0.2)) * 5)] + (1- alpha[(int)((currentValue / (Irange + 0.2)) * 5)])*newalpha;

		lookDirection.normalize();
		time++;

	} while (newalpha < 1.0||lookDirection.x<300||lookDirection.y<300||lookDirection.z<300);
	//printf("%f \n", currentColor.x);
	return currentColor;

}*/
unsigned char medium(float v, unsigned char l, unsigned char h)
{
	if (v > h) return h;
	else if (v < l) return l;
	else return (unsigned char)v;
}
RGBDATA trace( Vec3f &rayorig, Vec3f &raydir, float *volume) {
	Vec3f Q;
	Vec3f finalColor = { 0,0,0 };
	int time = 0;
	float alpha = 0.0f;
	float currentA = 0.0f;
	const int maxpossiblestep = (int)sqrt(xdim*xdim + ydim*ydim) / INTERVAL;
	
	for (int k = 0; k <= maxpossiblestep; k++) {
		
		Q = rayorig + (raydir*INTERVAL)*(float)k;
		
		if (alpha > 1.0) {
			break;
		}
		int index = IndexVect(min(round(Q.x), xdim - 1), min(round(Q.y), ydim-1),min(round(Q.z), zdim - 1));
		float currentValue = volume[index];
		Vec3f currentColor = {0,0,0};
		if (currentValue > 80)
		{
			currentColor = { 0, 0, 255 };
			currentA = 0.3f;
		}
		else if (currentValue > 70) {
			currentColor = {0, 255, 0};
			currentA = 0.3f;
		}
		else if (currentValue > 60) {
			currentColor = { 255, 0, 0 };
			currentA = 0.4f;
		}
		else if (currentValue > 50)
		{
			currentColor = { 0, 255, 255 };
			currentA = 0.3f;
		}
		else
		{
			currentColor = { 255, 255, 255 };
			currentA = 0.2f;
		}

		alpha = alpha + (1 - alpha) * currentA;
		finalColor = finalColor + currentColor*currentA;
		time++;
		//finalColor = (finalColor) / (65535 - 0) * 255;
		//finalColor = finalColor / time;
	}
	finalColor = finalColor / time*1.6;
	return{ medium(finalColor.x, 0, 255), medium(finalColor.y, 0, 255), medium(finalColor.z, 0, 255) };
}

void render(RGBDATA* output, float *volume, unsigned outXdim, unsigned outYdim) {
	//printf("I'm here in render");
	matrix xAxis = { { 1.0f, 0.0f, 0.0f },
	{ 0.0f, cosf(XR), sinf(XR) },
	{ 0.0f, -sinf(XR), cosf(XR) }
	};

	matrix yAxis = { { cosf(YR), 0.0f, sinf(YR) },
	{ 0.0f,  1.0f, 0.0f },
	{ -sinf(YR), 0.0f, cosf(YR) }
	};

	matrix zAxis = { { cosf(ZR), sinf(ZR), 0.0f },
	{ -sinf(ZR), cosf(ZR), 0.0f },
	{ 0.0f, 0.0f, 1.0f }
	};

	matrix R = xAxis*yAxis*zAxis;
	Vec3f B = {0.0,0.0,0.0};
	Vec3f S0 = {0.0,0.0,d};
	Vec3f u0 = { 1.0,0.0,0.0 };
	Vec3f v0 = {0.0,1.0,0.0};
	Vec3f g = {0.0,0.0,1.0};
	Vec3f S = B - (g*d);
	Vec3f u = mul(u0,R);
	Vec3f v = mul(v0, R);
	Vec3f E = S - (u *(outXdim / 2.0))  - (v* (outYdim / 2.0)) ;
	Vec3f T = { xdim / 2.0, ydim / 2.0, zdim / 2.0};
	Vec3f P;
	g = mul(g,R);

	/*Vec3f color[5] = {Vec3f(1,1,1), Vec3f(1,0,0), Vec3f(0,1,0) ,Vec3f(0,0,1) ,Vec3f(0,0,0) };*/
	/*float invWidth = 1 / float(outXdim), invHeight = 1 / float(outYdim);
	float fov = 90, aspectratio = outXdim / float(outYdim);
	float angle = tan(M_PI * 0.5 * fov / 180.);*/
	/*Vec3f rayorig;
	Vec3f raydir;*/
	for (unsigned y = 0; y < outYdim; ++y) {
		for (unsigned x = 0; x < outXdim; ++x) {
			P = E + u*(float)x + v*(float)y;
			P = P + T;
			output[x + outXdim * y] = trace(P,g,volume);
		}
	}

	/*std::ofstream ofs("./heart.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << outXdim << " " << outYdim << "\n255\n";
	for (unsigned i = 0; i < outXdim * outYdim; ++i) {
		ofs << (unsigned char)(min(255, image[i].x)) <<
			(unsigned char)(min(255, image[i].y) ) <<
			(unsigned char)(min(255, image[i].z) );
	}
	ofs.close();
	delete[] image;*/
}


int main(int argc, char *argv[])
{
  //int xdim,ydim,zdim;
  float *volume;
  FILE *fp;
  float minraw;
  float maxraw;
  

  if (argc != 3){
    printf("Usage: test <input_filename> <output> \n");
    printf("       <input_filename>:   RAWIV file \n");
    printf("       <output>:   PPM file \n");
    exit(0);              
  }

  
  //printf("begin reading rawiv.... \n");
  ReadRawiv(&xdim,&ydim,&zdim,&volume,argv[1]);
  float min = 100.0;
  float max = 0.0;


  //printf("I'm here in main");
  printf("%d\n", xdim);
  printf("%d\n", sizeof(ydim));
  printf("%d\n", sizeof(zdim));
  //define you output dimension
  int outXdim = 128;
  int outYdim = 128;
  RGBDATA *output;
  output = (RGBDATA*)malloc(sizeof(RGBDATA)*outXdim*outYdim);
  
  for (int i = 0; i < 128 * 128; i++) {
	  cout << volume[i] << endl;
  }
 // render(output, volume, outXdim, outYdim);
  // Volume rendering here ;
  render(output, volume, outXdim, outYdim);

  /* Begin writing PPM.... This is your output picture */
  // You may need to scale your output image to [0, 255]
  // before you write it to the PPM format.

  printf("Begin writing PPM.... \n");
  if ((fp=fopen(argv[2], "wb")) == NULL){
     printf("write ppm error....\n");
     exit(0);
   }
  WritePPM(outXdim, outYdim, output, fp);
  //printf("everything finish!");
  free(output);
  free(volume);
  //printf("everything finish!");
  //system("pause");
  return(0);
}

