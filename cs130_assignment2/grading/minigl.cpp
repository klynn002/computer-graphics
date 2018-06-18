/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h
//added
vec3 color;
MGLpoly_mode curtype;
mat4 projection;
MGLmatrix_mode cur_mat_mode;
mat4 modelview;
mat4 Id = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
struct Vertex{
	vec4 position;
	vec3 vercolor;
};
struct Triangle{
	Vertex A;
	Vertex B;
	Vertex C;
};
std::vector<Vertex> verlist;
std::vector<Triangle> trilist;
std::vector<mat4> proj_stack = {Id};
std::vector<mat4> model_stack = {Id};

MGLfloat areaof(vec2 a, vec2 b, vec2 c){
  float result = a[0]*(b[1] - c[1]);
  result += a[1]*(c[0]- b[0]);
  result += (b[0]*c[1] - b[1]*c[0]);
  return result;
}
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}
//added
void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data){
	
  
    float tri_area, alpha, beta, gamma; 
    
    MGLfloat fia = (tri.A.position[0]/tri.A.position[3] + 1)*width/2-0.5; 
    MGLfloat fja = (tri.A.position[1]/tri.A.position[3] + 1)*height/2-0.5; 
    MGLfloat fib = (tri.B.position[0]/tri.B.position[3] + 1)*width/2-0.5; 
    MGLfloat fjb = (tri.B.position[1]/tri.B.position[3] + 1)*height/2-0.5; 
    MGLfloat fic = (tri.C.position[0]/tri.C.position[3] + 1)*width/2-0.5; 
    MGLfloat fjc = (tri.C.position[1]/tri.C.position[3] + 1)*height/2-0.5;
    vec2 a(fia, fja);
    vec2 b(fib, fjb);   
    vec2 c(fic, fjc);     
    
    tri_area = areaof(a,b,c); 
 
    for(int i = 0; i < width; ++i) { 
       for (int j = 0; j < height; ++j) {
           vec2 p; 
           p[0] = i; 
           p[1] = j; 
   
           alpha = areaof(p,b,c)/tri_area; 
           beta = areaof(a,p,c)/tri_area; 
           gamma = areaof(a,b,p)/tri_area; 

           if ((alpha >= 0 && beta >=0 && gamma >= 0) 
              && (gamma >= 0 && gamma <= 1)){
              data[i+j*width] = Make_Pixel(tri.A.vercolor[0]*255, tri.A.vercolor[1]*255,tri.A.vercolor[2]*255); 
           }
       }
    }
}
mat4& current_matrix(){
    if(cur_mat_mode == MGL_PROJECTION){
      return projection;
    }
    else{
      return modelview;
    }
}
mat4& top_of_active_matrix_stack(){
  if(cur_mat_mode == MGL_PROJECTION){
    if(!proj_stack.empty()){
      return proj_stack.back();
    }
  }
  else{
    if(!model_stack.empty()){
      return model_stack.back();
    }
  }
  return current_matrix();
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
  Make_Pixel(0,0,0);
  for(unsigned int i = 0; i < trilist.size(); ++i){
      Rasterize_Triangle(trilist.at(i), width, height, data);
  }
  trilist.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	curtype = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(curtype == MGL_TRIANGLES){
		int vernum = (verlist.size() / 3) * 3;
		for(int i = 0; i < vernum; ++i){
			if((i + 1) % 3 == 0){
				Triangle t;
				t.A = verlist.at(i-2);
				t.B = verlist.at(i-1);
				t.C = verlist.at(i);
				trilist.push_back(t);
			}
		}
	}
	else if(curtype == MGL_QUADS){
		int vernum = (verlist.size() / 4) * 4;
		for(int i = 0; i < vernum; ++i){
			if((i + 1) % 4 == 0){
				Triangle t1;
				Triangle t2;
				t1.A = verlist.at(i-3);
				t1.B = verlist.at(i-2);
				t1.C = verlist.at(i-1);
				trilist.push_back(t1);
				t2.A = verlist.at(i-3);
				t2.B = verlist.at(i-1);
				t2.C = verlist.at(i);
				trilist.push_back(t2);
			}
		}
	}
	verlist.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	Vertex ver3;
  vec4 p = vec4(x,y,z,1);
  p = proj_stack.back() * model_stack.back() * p;
	ver3.position = p;
	ver3.vercolor = color;
  
	
	verlist.push_back(ver3);
	
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  cur_mat_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
  if(cur_mat_mode == MGL_PROJECTION){
    if(!proj_stack.empty()){
      proj_stack.push_back(proj_stack.back());
    }
  }
  else{
    if(!model_stack.empty()){
      model_stack.push_back(model_stack.back());
    }
  }
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{ 
  if(cur_mat_mode == MGL_PROJECTION){
    if(proj_stack.empty()){
      proj_stack.pop_back();
    }
  }
  else{
    if(!model_stack.empty()){
      model_stack.pop_back();
    }
  }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  mat4 id_matrix = {{1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1}};
  mat4& curr_stack = top_of_active_matrix_stack();
  curr_stack = id_matrix;
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  mat4 trans_matrix = {{1,0,0,0,
                        0,1,0,0,
                        0,0,1,0,
                        x,y,z,1}}; 
  mat4& curr_stack = top_of_active_matrix_stack();
  curr_stack = trans_matrix * curr_stack;
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
  float c = cos(angle*M_PI/180);
  float s = sin(angle*M_PI/180);
  vec3 v = vec3(x,y,z).normalized();
  x = v[0];
  y = v[1];
  z = v[2];
  
  mat4 rotate_matrix = {{x*x*(1-c)+c,x*y*(1-c)-z*s, x*z*(1-c)+y*s, 0,
			                   y*x*(1-c)+z*s,y*y*(1-c)+c, y*z*(1-c)-x*s, 0,
   			                 x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c, 0,
			                   0,0,0,1}};
  mat4& curr_stack = top_of_active_matrix_stack();
  curr_stack = rotate_matrix * curr_stack;
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  mat4 scale_matrix = {{ x, 0, 0, 0, 
                         0, y, 0, 0, 
                         0, 0, z, 0,
		                     0, 0, 0, 1}};
  mat4& curr_matrix =top_of_active_matrix_stack();
  curr_matrix = scale_matrix * curr_matrix;
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
  float A,B,C,D;
  A = (right+left)/(right-left);
  B = (top+bottom)/(top-bottom);
  C = -(far+near)/(far-near);
  D = -2*(far*near)/(far-near);
  
  mat4 frus_matrix = {{(2*near)/(right-left),0,0,0, 
                        0,(2*near)/(top-bottom),0,0,
                        A,B,C,-1, 
                        0, 0, D, 0}}; 
  mat4& curr_matrix = top_of_active_matrix_stack();
  curr_matrix = frus_matrix * curr_matrix;
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    mat4 ortho_matrix = {{2/(right - left),0,0,0,
                          0, 2/(top - bottom),0,0,
                          0,0,-2/(far-near),0,
                          -((right + left)/(right - left)), -((top + bottom)/(top - bottom)), 
                          -((far + near)/(far- near)), 1}};
    mat4& curr_matrix = top_of_active_matrix_stack();
    curr_matrix = ortho_matrix * curr_matrix;
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color = vec3(red, green, blue);
}
