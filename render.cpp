// Victor Priano
// Computer Graphics (Phong Model Shading)

////////////////////////////////////////////////////////////////////////
// A simple wrapper for to store 3D vectors
// This file is not my original work but was supplied to us.
// Functionality: This struct allows me to store the coordinates
//                of any vertex specified on any subdivision level
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;

struct Vector3
{
        float x,y,z;

        Vector3() : x(0.0), y(0.0), z(0.0)
        {}

        Vector3(float x, float y, float z)
                : x(x), y(y), z(z)
        {}

        Vector3(const Vector3 & v)
                : x(v.x), y(v.y), z(v.z)
        {}
        
        Vector3 operator+(const Vector3 & rhs) const
        { return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
        Vector3 operator-(const Vector3 & rhs) const
        { return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
        Vector3 operator*(float rhs) const
        { return Vector3(x * rhs, y * rhs, z * rhs); }
        Vector3 operator/(float rhs) const
        { return Vector3(x / rhs, y / rhs, z / rhs); }
        Vector3 operator+=(const Vector3 & rhs)
        { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
        Vector3 operator-=(const Vector3 & rhs)
        { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
        Vector3 operator*=(float rhs)
        { x *= rhs; y *= rhs; z *= rhs; return *this; }
        Vector3 operator/=(float rhs)
        { x /= rhs; y /= rhs; z /= rhs; return *this; }

        float magnitude() const
        { return sqrt(x * x + y * y + z * z); }
        void normalize()
        { *this /= magnitude(); }
        float dot(const Vector3 & rhs) const
        {
                return x * rhs.x + y * rhs.y + z * rhs.z;
        }
        Vector3 cross(const Vector3 & rhs) const
        {
                return Vector3(y * rhs.z - z * rhs.y,
                                        z * rhs.x - x * rhs.z,
                                        x * rhs.y - y * rhs.x);
        }
        
        void print() 
        { cout << x << " " << y << " " << z << endl, cout.flush(); }
};
////////////////////////////////////////////////////////////////////////
// 						           phong shading main                 		      //
////////////////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <limits>

// define constants to help with calculations
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
#define WEIGHT_1 3.0/8.0
#define WEIGHT_2 1.0/8.0
#define PI 3.14159265
#define INF numeric_limits<float>::infinity()

// variables dealing directly with the triangle mesh data
// vector structs to hold the coordinates and connectivity amongst them
vector< Vector3 > coordinates;
vector< Vector3 > vertices;
vector< vector<int> > connect;
vector< Vector3 > new_coordinates;
vector< Vector3 > new_vertices;
vector< Vector3 > tmp_coordinates;
vector< Vector3 > tmp_vertices;
vector< Vector3 > update_coordinates;
vector< Vector3 > update_vertices;
vector< Vector3 > normals;
int num_coordinates = 0;
int faces = 0;

// since OpenGl runs an infinite loop, I need flags to activate 
// certain functions
int view=0;
int divide=4;
int test_print=0;
int once=0;
int old_vertices=0;
int subdivision_level=1;
int flag_translate=0;
int flag_rotate=0;
int flag_render=0;
Vector3 center;

////////////////////////////////////////////////////////////////////////
// A storage container for an edge that specifies what vertex is
// adjacent and opposite of the edge
// Functionality: easier access to vertices that give weight to a
//                certain edge
////////////////////////////////////////////////////////////////////////
struct edges 
{
	int opp1;
	int opp2;
	int adj1;
	int adj2;
};
vector <edges> edgez;
////////////////////////////////////////////////////////////////////////
// Renders a quad at cell (x, y) with dimensions CELL_LENGTH          //
////////////////////////////////////////////////////////////////////////
void renderPixel(int x, int y, int z, float r = 0.0, float g = 0.0, float b = 0.0)
{
	// ...
	// Complete this function
	// ...
	glBegin(GL_POINTS);
	glColor3f(r,g,b);
	glVertex3i(x,y,0);	
	glEnd();
}
////////////////////////////////////////////////////////////////////////
//                      debugging                                     //
////////////////////////////////////////////////////////////////////////
// Simple debugging helper functions
// Functionality: Renders a point at (x,y) in yellow to
//                distinguish the newly inserted coordinate
void test_point(int x, int y)
{
	glPointSize(4.0);
	glColor3f(1.0, 1.0, 0.0);
	glEnable(GL_POINT_SMOOTH);
	renderPixel(x, y, 0);
	glPointSize(1.0);
	glColor3f(1.0, 1.0, 1.0);
	glDisable(GL_POINT_SMOOTH);
}
// Functionality: Prints to the terminal the values of the vectors
//                containing the coordinates, face list, and adjcency
//				  list, mostly used for debugging purposes
void print()
{
	cout << "VERTICES\n", cout.flush();
	for (int i=0; i<new_vertices.size();++i)
		cout << new_vertices[i].x << " " 
		     << new_vertices[i].y << " " 
		     << new_vertices[i].z << endl, cout.flush();
	cout << "COORDINATES\n", cout.flush();
	for (int j=0; j<new_coordinates.size();++j)
		cout << "( " << new_coordinates[j].x << " " 
		     << new_coordinates[j].y << " " << new_coordinates[j].z 
		     << " )" << endl, cout.flush();
	cout << "ADJACENCY LIST\n", cout.flush();
	for (int k=0; k<connect.size();++k)
	{	
		cout << connect[k][0] << ": ", cout.flush();
		for (int l=1; l<connect[k].size();++l)
			cout << connect[k][l] << " ", cout.flush();
		cout << endl, cout.flush();
	}
}
// Functionality : Calculates the centroid of the mesh which
//                 will help determining the right direction of the
//                 normal
float centroid(vector<Vector3> c, int n);
Vector3 find_centroid(vector<Vector3>c)
{
	float x_=centroid(c,0);
	float y_=centroid(c,1);
	float z_=centroid(c,2);
	renderPixel(x_*8,y_*8,0,1.0,0.0,0.0);
	return Vector3(x_,y_,z_);
}
// Functionality: calculates the smallest x coordinate in the mesh
float xmin(vector<Vector3>c)
{
	float min=c[0].x;
	for (int i=1; i<c.size();++i) { if (c[i].x<min) { min=c[i].x; }}
	return min;
}
// Functionality: calculates the largest x coordinate in the mesh
float xmax(vector<Vector3>c)
{
	float max=c[0].x;
	for (int i=1; i<c.size();++i) { if (c[i].x>max) { max=c[i].x; }}
	return max;
}
// Functionality: calculates the smallest y coordinate in the mesh
float ymin(vector<Vector3>c)
{
	float min=c[0].y;
	for (int i=1; i<c.size();++i) { if (c[i].y<min) { min=c[i].y; }}
	return min;
}
// Functionality: calculates the largest y coordinate in the mesh
float ymax(vector<Vector3>c)
{
	float max=c[0].y;
	for (int i=1; i<c.size();++i) { if (c[i].y>max) { max=c[i].y; }}
	return max;
}
// Functionality: calculates the smallest z coordinate in the mesh
float zmin(vector<Vector3>c)
{
	float min=c[0].z;
	for (int i=1; i<c.size();++i) { if (c[i].z<min) { min=c[i].z; }}
	return min;
}
// Functionality: calculates the largest z coordinate in the mesh
float zmax(vector<Vector3>c)
{
	float max=c[0].z;
	for (int i=1; i<c.size();++i) { if (c[i].z>max) { max=c[i].z; }}
	return max;
}
////////////////////////////////////////////////////////////////////////
//                    Graphics rendering                              //
////////////////////////////////////////////////////////////////////////
// Rendering algorithms (completed in lab)
// Functionality: DDA renders a line from (x1, y1) to (x2, y2)
void DDA (int x1, int y1, int x2, int y2)
{
	int dx = x2-x1;
	int dy = y2-y1;
	float m = (float)dy/(float)dx;
	int steps = (abs(dx) > abs(dy)) ? abs(dx) : abs(dy);			
	float x=x1, y=y1;

	float xinc = dx/(float)steps;
	float yinc = dy/(float)steps;
	
	for (int i = 0; i < steps; ++i)
	{
		x+=xinc;
		y+=yinc;
		renderPixel(x, y, 0);
	}
}
// Functionality: Renders part of a circle using an initial start (x,y)
// 				  with radius r and moving along the arc calculated
//			      through a formula
void MidpointCircle(int x, int y, int r)
{
	int i=0, j=r;
	float midpoint = 0;
	while (i <= j)
	{
		renderPixel(x+i, y+j, 0);
		renderPixel(x-i, y+j, 0);
		renderPixel(x+i, y-j, 0);
		renderPixel(x-i, y-j, 0);
		renderPixel(x+j, y+i, 0);
		renderPixel(x-j, y+i, 0);
		renderPixel(x+j, y-i, 0);
		renderPixel(x-j, y-i, 0);
		++i;

		midpoint = (r*r) - (i*i) - ((j+0.5)*(j+0.5));
		j = (midpoint < 0) ? j-1 : j;
	}
}
////////////////////////////////////////////////////////////////////////
//                          read data                                 //
////////////////////////////////////////////////////////////////////////
// Reads a text file through command line
// Functionality: Parses the specified file and stores the files 
//                contents in the appropriate container
void read_file(char* filename)
{
	float data = 0.0;
	ifstream file;
	file.open(filename);
	if (file == NULL) 
	{ 
		cout << "read_file() failed: FILE DOES NOT EXIST.\n";
		exit(1); 
	}
	file >> data;
	num_coordinates = data;
	file >> data;
	faces = data;
	Vector3 point;
	for (int i=0; i < num_coordinates; ++i)
	{
		file >> data;
		point.x = data;
		file >> data;
		point.y = data;
		file >> data;
		point.z = data;
		coordinates.push_back(point);
	}
	Vector3 edge;
	for (int i=0; i < faces; ++i)
	{
		file >> data;
		edge.x = data;
		file >> data;
		edge.y = data;
		file >> data;
		edge.z = data;
		vertices.push_back(edge);
	}
}
////////////////////////////////////////////////////////////////////////
//                       helper functions                             //
////////////////////////////////////////////////////////////////////////
// Simple helper functions relating to the subdividing functions
// The applied weight to an existing vertex
//Function
float alpha(int n)
{
	float a = 1.0/4.0;
	float b = cos(2*PI/n)*a;
	b += WEIGHT_1;
	b *= b;
	float c = 5.0/8.0;
	c -= b;
	return (1.0/float(n))*c;
}

float beta(int n)
{
	return 1.0 - (float)n*alpha(n);
}
////////////////////////////////////////////////////////////////////////
//                     helper functions                               //
////////////////////////////////////////////////////////////////////////
bool is_contained_int(vector<int> v, float n)
{
	for (int i=0; i< v.size(); ++i) if ((float)v[i] == n) return true;
	return false;
}

bool is_contained_Vector3(vector<Vector3> v, Vector3 t)
{
	for (int i=0; i< v.size(); ++i) if (v[i].x == t.x && v[i].y == t.y && v[i].z == t.z) return true;
	return false;
}

int find(vector<Vector3>v, Vector3 t)
{
	for (int i=0; i< v.size(); ++i) if (v[i].x == t.x && v[i].y == t.y && v[i].z == t.z) return i;
	return -1;
}
////////////////////////////////////////////////////////////////////////
//                      main subdivide functions                      //
////////////////////////////////////////////////////////////////////////
// Creates the adjacency list
// Functionality: Updates the connect vector which stores a vertex 
//			 	  index and the indices of all vertices that connect to
//                that specific vertex
void connectivity(vector<Vector3> v, vector<Vector3> c)
{
	vector<int> tmp;
	connect.clear();
	for(int i=0; i < c.size(); ++i) tmp.push_back(i), connect.push_back(tmp), tmp.clear();
	for (int j=0; j < c.size(); ++j)
	{
		for (int k=0; k < v.size(); ++k)
		{
			if (v[k].x == connect[j][0])
			{
				if (!is_contained_int(connect[j], v[k].y)) connect[j].push_back(v[k].y);
				if (!is_contained_int(connect[j], v[k].z)) connect[j].push_back(v[k].z);
			}
			if (v[k].y == connect[j][0])
			{
				if (!is_contained_int(connect[j], v[k].x)) connect[j].push_back(v[k].x);
				if (!is_contained_int(connect[j], v[k].z)) connect[j].push_back(v[k].z);
			}
			if (v[k].z == connect[j][0])
			{
				if (!is_contained_int(connect[j], v[k].x)) connect[j].push_back(v[k].x);
				if (!is_contained_int(connect[j], v[k].y)) connect[j].push_back(v[k].y);
			}
		}
	}
}
// Functionality: applies weight to existing vertices using neighbors
void update_even(vector<Vector3> c)
{
	float alpha_constant=0.0;
	float beta_constant=0.0;
	for (int i=0; i < old_vertices; ++i)
	{
		Vector3	update = c[i];
		Vector3 tmp;
		alpha_constant = alpha(connect[i].size()-1);
		beta_constant = beta(connect[i].size()-1);
		update*= beta_constant;
		for (int j=1; j < connect[i].size(); ++j) { tmp+=c[connect[i][j]]; }
		tmp*=alpha_constant;
		update+=tmp;
		new_coordinates[i] = update;
	}
}
// Functionality: returns the midpoint that has been already weighed
Vector3 update_odd(vector<Vector3>c, int n, int x)
{
	Vector3 odd;
	for (int i=0; i<edgez.size();++i)
	{
		if ( (edgez[i].adj1==n && edgez[i].adj2==x) || (edgez[i].adj1==x && edgez[i].adj2==n) )
		{
			odd.x=WEIGHT_1*c[n].x+WEIGHT_1*c[x].x+WEIGHT_2*c[edgez[i].opp1].x+WEIGHT_2*c[edgez[i].opp2].x;
			odd.y=WEIGHT_1*c[n].y+WEIGHT_1*c[x].y+WEIGHT_2*c[edgez[i].opp1].y+WEIGHT_2*c[edgez[i].opp2].y;
			odd.z=WEIGHT_1*c[n].z+WEIGHT_1*c[x].z+WEIGHT_2*c[edgez[i].opp1].z+WEIGHT_2*c[edgez[i].opp2].z;
			return odd;
		}
	}
	return odd;
}
// Functionality: groups edges into a premad edge struct that stores
//                each edges adjacent and opposite points
void group_edges(vector<Vector3> v)
{
	float side_1=0.0, side_2=0.0,side_3=0.0;
	float next1=0.0, next2=0.0, next3=0.0;
	edgez.clear();
	for (int i=0; i < v.size()-1; ++i)
	{
		side_1 = v[i].x;
		side_2 = v[i].y;
		side_3 = v[i].z;
		for (int j=i+1; j < v.size(); ++j)
		{
			edges line;
			next1 = v[j].x;
			next2 = v[j].y;
			next3 = v[j].z;
			if ( (side_1 == next1 || side_1 == next2 || side_1 == next3) &&
			     (side_2 == next1 || side_2 == next2 || side_2 == next3) )
			{
				line.adj1 = side_1, line.adj2 = side_2, line.opp1 = side_3;
				if (next1 != side_1 && next1 != side_2) line.opp2 = next1;
				else if (next2 != side_1 && next2 != side_2) line.opp2 = next2;
				else line.opp2 = next3;
				edgez.push_back(line);
	 	     }
			else if ( (side_1 == next1 || side_1 == next2 || side_1 == next3) && 
			     	  (side_3 == next1 || side_3 == next2 || side_3 == next3))
			{
				line.adj1 = side_1, line.adj2 = side_3, line.opp1 = side_2;
				if (next1 != side_1 && next1 != side_3) line.opp2 = next1;
				else if (next2 != side_1 && next2 != side_3) line.opp2 = next2;
				else line.opp2 = next3;
				edgez.push_back(line);
			}
			else if ( (side_2 == next1 || side_2 == next2 || side_2 == next3) && 
				      (side_3 == next1 || side_3 == next2 || side_3 == next3))
			{
				line.adj1 = side_2, line.adj2 = side_3, line.opp1 = side_1;
				if (next1 != side_2 && next1 != side_3) line.opp2 = next1;
				else if (next2 != side_2 && next2 != side_3) line.opp2 = next2;
				else line.opp2 = next3;
				edgez.push_back(line);
			}
			else {}
		}
	}
}
// Functionality: The meat of the loop subdivision code. It generates
//                the midpoint location and ask for its new position.
//                It updates the coordinates and vertices vectors
void subdivide(vector<Vector3> v, vector<Vector3> c)
{
	new_vertices.clear();
	new_coordinates.clear();
	new_coordinates=c;
	old_vertices=c.size();
	Vector3 update, mid1, mid2, mid3;
	for (int i=0; i<v.size();++i)
	{

		mid1 = update_odd(c, v[i].x, v[i].y);
		if (!is_contained_Vector3(new_coordinates, mid1)) new_coordinates.push_back(mid1);
		mid2 = update_odd(c, v[i].x, v[i].z);
		if (!is_contained_Vector3(new_coordinates, mid2)) new_coordinates.push_back(mid2);
		mid3 = update_odd(c, v[i].y, v[i].z);
		if (!is_contained_Vector3(new_coordinates, mid3)) new_coordinates.push_back(mid3);
		update.x=v[i].x;
		update.y=find(new_coordinates, mid1);
		update.z=find(new_coordinates, mid2);
		new_vertices.push_back(update);
		update.x=v[i].y;
		update.y=find(new_coordinates, mid1);
		update.z=find(new_coordinates, mid3);
		new_vertices.push_back(update);
		update.x=v[i].z;
		update.y=find(new_coordinates, mid2);
		update.z=find(new_coordinates, mid3);
		new_vertices.push_back(update);
		update.x=find(new_coordinates, mid1);
		update.y=find(new_coordinates, mid2);
		update.z=find(new_coordinates, mid3);
		new_vertices.push_back(update);
	}
}
////////////////////////////////////////////////////////////////////////
// 				     orientations of image			                  //
////////////////////////////////////////////////////////////////////////
// Functionality: calculates the average of the coordinates and where
//                the center of mass is
float centroid (const vector<Vector3>c, int n)
{
	float xmiddle=0.0; float ymiddle=0.0; float zmiddle =0.0;
	for(int i=0; i<c.size();++i) { xmiddle+=c[i].x, ymiddle+=c[i].y, zmiddle+=c[i].z; }
	return (n==0) ? xmiddle/(float)c.size():(n==1)?ymiddle/(float)c.size():zmiddle/(float)c.size();
}
// Functionality: this function translates the mesh either
//                horizontally or vertically on a predefined key press
void translate(vector<Vector3>&c, int t_orientation, int t_direction, int v)
{
	Vector3 tmp;
	for (int i=0; i<c.size();++i)
	{
		tmp.x=c[i].x, tmp.y=c[i].y, tmp.z=c[i].z;
		if (v==0)
		{
			if (t_orientation) { c[i].x = (t_direction) ?  tmp.x+1.0 : tmp.x-1.0; }
			else { c[i].y = (t_direction) ?  tmp.y+1.0 : tmp.y-1.0; }
		}
		if (v==1)
		{
			if (t_orientation) { c[i].x = (t_direction) ?  tmp.x+1.0 : tmp.x-1.0; }
			else { c[i].z = (t_direction) ?  tmp.z+1.0 : tmp.z-1.0; }
		}
		if (v==2)
		{
			if (t_orientation) { c[i].y = (t_direction) ?  tmp.y+1.0 : tmp.y-1.0; }
			else { c[i].z = (t_direction) ?  tmp.z+1.0 : tmp.z-1.0; }
		}
	}
}
// Functionality: the Quaternion struct helps with rotations
struct Quaternion
{
	float w,x,y,z;
	Quaternion() : w(0.0), x(0.0), y(0.0), z(0.0) {}
	Quaternion(float angle, float _x, float _y, float _z) 
	{ w=cosf(angle/2), x=_x*sinf(angle/2),y=_y*sinf(angle/2),z=_z*sinf(angle/2); }
};
// rotating global helper constants
float degree=5.0;
float radians = degree*0.0174532925;
const Quaternion pos_x(radians,1,0,0);
const Quaternion neg_x(radians,-1,0,0);
const Quaternion pos_y(radians,0,1,0);
const Quaternion neg_y(radians,0,-1,0);
const Quaternion pos_z(radians,0,0,1);
const Quaternion neg_z(radians,0,0,-1);
// Functionality: sets up the quaternion matrix based on an angle, w 
//                and three coordinates of an axis, x y and z
vector<Vector3> rotation_matrix(Quaternion q)
{
	vector<Vector3>matrix;
	Vector3 a,b,c;
	a.x=1-2*q.y*q.y-2*q.z*q.z;
	a.y=2*q.x*q.y-2*q.w*q.z;
	a.z=2*q.x*q.z+2*q.w*q.y;
	b.x=2*q.x*q.y+2*q.w*q.z;
	b.y=1-2*q.x*q.x-2*q.z*q.z;
	b.z=2*q.y*q.z+2*q.w*q.x;
	c.x=2*q.x*q.z-2*q.w*q.y;
	c.y=2*q.y*q.z-2*q.w*q.x;
	c.z=1-2*q.x*q.x-2*q.y*q.y;
	
	matrix.push_back(a), matrix.push_back(b), matrix.push_back(c);
	return matrix;
}
// Functionality: this function rotates the object along a specific
//                axis based on which key press the user entered
void rotate(vector<Vector3>&c, int r_direction, int r_orientation, int v)
{
	vector<Vector3> rotation;
	Vector3 tmp;
	int axis=0;
	Quaternion q;
	if (v==0)
	{
			if (r_orientation) { q = (r_direction) ?  pos_x : neg_x, axis=0; }
			else { q = (r_direction) ?  pos_y : neg_y, axis=1; }
	}
	if (v==1)
	{
			if (r_orientation) { q = (r_direction) ?  pos_x : neg_x, axis=0; }
			else { q = (r_direction) ?  pos_z : neg_z, axis=2; }
	}
	if (v==2)
	{
			if (r_orientation) { q = (r_direction) ?  pos_y : neg_y, axis=1; }
			else { q = (r_direction) ?  pos_z : neg_z, axis=2; }
	}
	
	rotation=rotation_matrix(q);
	float average=centroid(c,axis);
	
	for(int i=0; i<c.size();++i)
	{
		tmp=c[i];
		tmp.x-=average, tmp.y-=average, tmp.z-=average;
		c[i].x=rotation[0].dot(tmp);
		c[i].y=rotation[1].dot(tmp);
		c[i].z=rotation[2].dot(tmp);
		c[i].x+=average, c[i].y+=average, c[i].z+=average;
	}
	
}
////////////////////////////////////////////////////////////////////////
//                      lighting/shading                              //
////////////////////////////////////////////////////////////////////////
// Phoung lighting global constants which are predefined and assumed
const Vector3 LIGHTSOURCE(800.0,800.0,0.0);
const Vector3 VIEWSOURCE (0.0,800.0,0.0);
const Vector3 material(0.0,0.0,1.0);
const Vector3 ambient(0.4,0.4,0.4);
const Vector3 diffuse(0.7,0.7,0.7);
const Vector3 specular(1.0,1.0,1.0);
vector< vector<int> >shared;
// Functionality: Groups the indicies that share a specific vertex.
//                The purpose is to find the averaged normal for that
//                vertex based on its neighbor normals
void group_vertex(vector<Vector3>v)
{
	vector<int> tmp;
	shared.clear();
	for(int i=0; i < v.size(); ++i) tmp.push_back(i), shared.push_back(tmp), tmp.clear();
	for (int i=0; i<v.size();++i)
	{
		for (int j=0;j<v.size();++j) { if (v[j].x==i || v[j].y==i || v[j].z==i) { shared[i].push_back(j); } }
	}
}
// Functionality: calculates the pixel normal and the averaged normal
void calc_norms(vector<Vector3>c, vector<Vector3>v)
{
	normals.clear();
	Vector3 pt1,pt2,pt3,hold, U,V;
	group_vertex(v);
	int tmp=0;
	for(int i=0;i<v.size();++i)
	{
		Vector3 cp;
		for (int j=1; j<shared[i].size();++j)
		{
			tmp=shared[i][j];
			pt1=c[v[tmp].x],pt2=c[v[tmp].y],pt3=c[v[tmp].z];
			U=pt2-pt1;
			V=pt3-pt1;
			hold=U.cross(V);
			cp+=hold;
		}
		cp/=shared[i].size();
		cp.normalize();
		normals.push_back(cp);
	}
}
// Functionality: calculates the surface normal for one face. It 
//                linear interpolates among the triangles and averages
//                out the normals
float normal_beta=0.0, normal_gamma=0.0;
Vector3 surface_normal(vector<Vector3>v, int n)
{
	float normal_alpha=1.0-normal_beta-normal_gamma;
	return normals[v[n].x]*normal_alpha+normals[v[n].y]*normal_beta+normals[v[n].z]*normal_gamma;
}
// Functionality: finds the angle between two vectors. The return 
//                type is in degrees and not radians
float angle(Vector3 u, Vector3 v)
{
	float num=u.dot(v);
	float den=u.magnitude()*v.magnitude();
	return acosf(num/den)*57.2957795;
}
// Functionality: calculates the determinant of a 2 by 2 matrix
float det2x2(float a, float d, float b, float c)
{
	return a*d-b*c;
}
// Functionality: calculates the determinant of a 3 by 3 matrix
float det(float a[3][3] )
{
	return a[0][0]*det2x2(a[1][1],a[2][2],a[1][2],a[2][1])-
	       a[0][1]*det2x2(a[1][0],a[2][2],a[2][0],a[1][2])+
	       a[0][2]*det2x2(a[1][0],a[2][1],a[2][0],a[1][1]);
}
// Functionality: barycentric coordinates calculate whether the 
//                ray in the Rd direction hits any faces of the mesh
Vector3 Rdx = Vector3(1, 0, 0);
Vector3 Rdy = Vector3(0, 1, 0);
Vector3 Rdz = Vector3(0, 0, 1);
float intersect(Vector3 a, Vector3 b, Vector3 c, Vector3 Ro, Vector3 Rd)
{
	float return_t=0;
	Vector3 Eb = b-a;
	Vector3 Ec = c-a;
	float T[3][3];
	float A[3][3];
	float B[3][3];
	float C[3][3];

	T[0][0]= -Eb.x; T[0][1]=-Ec.x; T[0][2]=a.x-Ro.x;
	T[1][0]= -Eb.y; T[1][1]=-Ec.y; T[1][2]=a.y-Ro.y;
	T[2][0]= -Eb.z; T[2][1]=-Ec.z; T[2][2]=a.z-Ro.z;
	
	C[0][0]= -Eb.x; C[0][1]=a.x-Ro.x; C[0][2]=Rd.x;
	C[1][0]= -Eb.y; C[1][1]=a.y-Ro.y; C[1][2]=Rd.y;
	C[2][0]= -Eb.z; C[2][1]=a.z-Ro.z; C[2][2]=Rd.z;
	
	B[0][0]= a.x-Ro.x; B[0][1]=-Ec.x; B[0][2]=Rd.x;
	B[1][0]= a.y-Ro.y; B[1][1]=-Ec.y; B[1][2]=Rd.y;
	B[2][0]= a.z-Ro.z; B[2][1]=-Ec.z; B[2][2]=Rd.z;
	
	A[0][0]= -Eb.x; A[0][1]=-Ec.x; A[0][2]=Rd.x;
	A[1][0]= -Eb.y; A[1][1]=-Ec.y; A[1][2]=Rd.y;
	A[2][0]= -Eb.z; A[2][1]=-Ec.z; A[2][2]=Rd.z;
	if (det(A) != 0.0)
	{	
		return_t = det(T)/det(A);
		normal_beta = det(B)/det(A);
		normal_gamma = det(C)/det(A);
		return (return_t>0.0 && normal_beta>0.0 && normal_gamma>0.0 && normal_beta+normal_gamma<1.0) ? return_t : INF;
	}
	return INF;
}
// ray tracer global constants
float t=INF,shadow_t=INF, prev=INF;
int color=0;
float LN=0.0, VN=0.0;
// Functionality: traverses the 800 by 800 pixel window and updates
//                the color of the material based on the illumination
//                and material of the object
void ray_tracing(vector<Vector3>c, vector<Vector3>v, vector<Vector3> n, int V)
{
	int boundW1=0, boundW2=0, boundH1=0, boundH2=0;
	Vector3 RD, Ro, R;
	if (V==0)
	{
		boundW1=xmin(c)*8;
		boundW2=xmax(c)*8;
		boundH1=ymin(c)*8;
		boundH2=ymax(c)*8;
		RD=Rdz;
	}
	else if (V==1)
	{
		boundW1=xmin(c)*8;
		boundW2=xmax(c)*8;
		boundH1=zmin(c)*8;
		boundH2=zmax(c)*8;
		RD=Rdy;
	}
	else
	{
		boundW1=ymin(c)*8;
		boundW2=ymax(c)*8;
		boundH1=zmin(c)*8;
		boundH2=zmax(c)*8;
		RD=Rdx;
	}
	
	float sanity=0.0, shininess=0.0;
	Vector3 illumination, l, vw, middle, check, diff, spec;
	middle=find_centroid(c);
	for (int x=0; x<= WINDOW_WIDTH;++x)
	{
		for (int y=0;y<=WINDOW_HEIGHT;++y)
		{
			if (V==0) { Ro=Vector3(float(x),float(y), 0.0); }
			else if (V==1) { Ro=Vector3(float(x), 0.0, float(y)); }
			else { Ro=Vector3(0.0, float(x), float(y)); }
			Vector3 Fn;
			vw=RD*-1;
			diff=diffuse, diff.normalize();
			spec=specular, spec.normalize();
			color=3;
			if ( x<boundW1 || x>boundW2 || y<boundH1 || y>boundH2 )
			{
				l=LIGHTSOURCE*1;
				l-=Ro, l.normalize();
				l.z=1.0;//, l.normalize();
				shadow_t=INF, prev=INF, color=3;
				for(int i=0;i<v.size();++i)
				{
					shadow_t=intersect(c[v[i].x]*8,c[v[i].y]*8, c[v[i].z]*8, Ro, l);
					if (shadow_t<prev) 
					{ 
						prev=shadow_t;
						color=2;
						break; 
					}
				}
			}
			else
			{
				l=LIGHTSOURCE*-1;
				l-=Ro, l.normalize();
				vw+=l, vw.normalize();
				t=INF, prev=INF, color=3;
				for (int i=0; i<v.size();++i)
				{
					t=intersect(c[v[i].x]*8,c[v[i].y]*8, c[v[i].z]*8, Ro, RD);
					if (t<prev)
					{
						color=0;
						Fn=surface_normal(v,i);
						prev=t;
					}
				}
				      
				//check=middle-Ro;
				//sanity=check.dot(Fn);
				//if (sanity<0) { Fn*=-1; }
				
				// diffuse
				LN=l.dot(Fn);
				LN=(LN<0)?0:(LN>1)?1:LN;
				diff*=LN;
			
				// specular
				R=Fn*LN, R*=2;
				R=l-R;
				R.normalize();
				VN=vw.dot(R);
				VN=(VN<0)?0:(VN>1)?1:VN;
				shininess=pow(VN, 200);
				spec*=shininess;
				
				illumination=ambient+diff+spec;
				//llumination.print();
			}
			
			switch (color)
			{
				case 0:
					renderPixel(x,y,0,material.x*illumination.x,material.y*illumination.y,material.z*illumination.z);
					break;
				case 2:
					renderPixel(x,y,0,0.0,0.0,0.0);
					break;
				default:
					renderPixel(x,y,0,ambient.x,ambient.y,ambient.z);
					break;
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
//                      print                                         //
////////////////////////////////////////////////////////////////////////
// Functionality:  prints whatever is in the coordinates vector based on
//                 the adjacency list
void print_mesh(vector<Vector3> c, int n)
{
	for (int i = 0; i < connect.size();++i)
		for (int j=1; j < connect[i].size(); ++j)
		{
			if (!n) DDA(c[i].x*8, c[i].y*8, c[connect[i][j]].x*8, c[connect[i][j]].y*8);
			else if (n==1) DDA(c[i].x*8, c[i].z*8,c[connect[i][j]].x*8, c[connect[i][j]].z*8);
			else DDA(c[i].y*8, c[i].z*8, c[connect[i][j]].y*8, c[connect[i][j]].z*8);
		}
}
// basic GL_render() call
void GL_render()
{
	glClear(GL_COLOR_BUFFER_BIT);

	center=find_centroid(new_coordinates);
	
	if (!flag_render) { print_mesh(new_coordinates, view); }
	else  { calc_norms(new_coordinates, new_vertices), ray_tracing(new_coordinates, new_vertices, normals, view); }
	glutPostRedisplay();
	
	glutSwapBuffers();
}
// Functionality: updates program based on keyboard input
void KeyboardInput(unsigned char key, int x, int y)
{
	switch (key)
	{
		//    T    //
		case 84:
		case 116:
			flag_translate=1;
			flag_rotate=0;
			break;
		//    A    //
		case 65:
		case 97:
			if (flag_translate) { translate(new_coordinates,1,0, view), translate(update_coordinates,1,0, view); }
			if (flag_rotate) { rotate(new_coordinates, 0,0, view); }
			break;
		//    D    //
		case 68:
		case 100:
			if (flag_translate) { translate(new_coordinates,1,1, view), translate(update_coordinates,1,1, view); }
			if (flag_rotate) { rotate(new_coordinates, 1,0, view); }
			break;
		//    W    //
		case 87:
		case 119:
			if (flag_translate) { translate(new_coordinates,0,1, view), translate(update_coordinates,0,1, view); }
			if (flag_rotate) { rotate(new_coordinates, 1, 1, view); }
			break;
		//    S    //
		case 83:
		case 115:
			if (flag_translate) { translate(new_coordinates,0,0, view), translate(update_coordinates,0,0, view); }
			if (flag_rotate) { rotate(new_coordinates, 0,1, view); }
			break;
		//    V    //
		case 86:
		case 118:
			subdivision_level=1;
			++view;
			view = (view > 2) ? 0 : view;
			new_coordinates.clear(), tmp_coordinates.clear();
			new_vertices.clear(), tmp_vertices.clear();
			new_coordinates=coordinates, tmp_coordinates=coordinates;
			new_vertices=vertices, tmp_vertices=vertices;
			connectivity(new_vertices, new_coordinates);
			group_edges(new_vertices);
			break;
		//    L    //
		case 76:
		case 108:
			++subdivision_level;
			if (subdivision_level<5)
			{
				tmp_coordinates=new_coordinates;
				tmp_vertices=new_vertices;
				subdivide(tmp_vertices, tmp_coordinates);
				update_even(tmp_coordinates);
				connectivity(new_vertices, new_coordinates);			
				group_edges(new_vertices);
			}
			else
			{ 
				subdivision_level=1;
				new_coordinates=update_coordinates;
				new_vertices=update_vertices; 
				connectivity(new_vertices, new_coordinates);			
				group_edges(new_vertices);
			}
			break;
		//    R    //
		case 82:
		case 114:
			flag_rotate=1;
			flag_translate=0;
			break;
		//    I    //
		case 73:
		case 105:
			flag_render=(flag_render)?0:1;
			break;
		default:
			break;
	}
}
////////////////////////////////////////////////////////////////////////
//                      initialize original                           // 
////////////////////////////////////////////////////////////////////////
void subdivide_init()
{
	new_vertices=vertices, tmp_vertices=vertices, update_vertices=vertices;
	new_coordinates=coordinates, tmp_coordinates=coordinates, update_coordinates=coordinates;
	connectivity(vertices, coordinates);
	group_edges(vertices);
}
////////////////////////////////////////////////////////////////////////
//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	// ...
	// Complete this function
	// ...
	glutCreateWindow(“Victor Priano“);

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// For the purposes of this lab, this is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glutDisplayFunc(GL_render);
	glutKeyboardFunc(KeyboardInput);
}
////////////////////////////////////////////////////////////////////////
//                       main code                                    //
////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{	
	read_file(argv[1]);
	subdivide_init();
	GLInit(&argc, argv);
	glutMainLoop();

	return 0;
}

