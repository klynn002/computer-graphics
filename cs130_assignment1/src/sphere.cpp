#include "sphere.h"
#include "ray.h"


// Determine if the ray intersects with the sphere
bool Sphere::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
    vec3 p = ray.endpoint - this->center;
    double b = 2 * dot(ray.direction, p);
    double c = dot(p,p) - (this->radius * this->radius);
    double delta = (b * b) - (4 * c);
    if(delta < 0){
		return false;
	}
	if(delta > 0){
		double t1 = (-1 * b) - sqrt(delta);
   t1 = t1/2.0;
		double t2 = (-1 * b) + sqrt(delta);
   t2 = t2/2.0;
   Hit h1;
   Hit h2;
		if(t2 > 0){
      if(t1 > 0){
        h1.object = this;
        h1.t = t1;
        h1.ray_exiting = false;
        hits.push_back(h1);
        }
        else{
        h1.object = this;
        h1.t = 0;
        h1.ray_exiting = false;
        hits.push_back(h1);
        }
      h2.object = this;
      h2.t = t2;
      h2.ray_exiting = true;
      hits.push_back(h2);
      return true;
	}
   else{
     if(delta == 0){
       double t = (-b / 2);
       if(t > 0){
       Hit h3;
       h3.object = this;
       h3.t = t;
       h3.ray_exiting = true;
       hits.push_back(h3);
       }
       return true;
     }
     else{
       return false;
       }
 
 /*if(delta == 0){
		Hit h1;
		h1.object = this;
		h1.t = ;
		h1.ray_exiting = false;
		return true;
	}*/
 }
    }
    return false;
}

vec3 Sphere::Normal(const vec3& point) const
{
    vec3 normal;
    // TODO: set the normal
    normal = point - this->center;
    normal = normal.normalized();
    return normal;
}
