#include "plane.h"
#include "ray.h"
#include <cfloat>


// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
bool Plane::
Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
     double top = dot(this->normal, (this->x1 - ray.endpoint));
    double bot = dot(this->normal, ray.direction);
    if(bot) { //Check for a denominator of 0
	    double t1 = top/bot;
	      if(t1 > 0) { //If the Hit is in front of the ray
		      Hit h;
   	      h.object = this;
		      h.t=t1;
		      h.ray_exiting = false;
		      hits.push_back(h);
		      return true;
	      }
	    else { //If the hit was behind the ray, or the ray starts in the plane;
		      return false;
	    }
    }

    return false;
}

vec3 Plane::
Normal(const vec3& point) const
{
    return normal;
}
