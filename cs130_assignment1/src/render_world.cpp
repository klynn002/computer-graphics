#include <vector>
#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"


Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3),disable_fresnel_reflection(false),disable_fresnel_refraction(false)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find the closest object of intersection and return the object that was
// intersected.  Record the Hit structure in hit.  If no intersection occurred,
// return NULL.  Note that in the case of a Boolean, the object returned will be
// the Boolean, but the object stored in hit will be the underlying primitive.
// Any intersection with t<=small_t should be ignored.
Object* Render_World::Closest_Intersection(const Ray& ray,Hit& hit)
{
    Object* closestObject = 0;
    double min_t = 999999999999999.0;
    for(unsigned int i = 0; i < objects.size(); ++i){
		std::vector<Hit> h;
		if(objects.at(i)->Intersection(ray, h)){
		  for(unsigned int j = 0; j < h.size(); ++j){
			  if(h.at(j).t < min_t){
           min_t = h.at(j).t;
           hit = h.at(j);
           closestObject = objects.at(i);              
       }
		  }
   }
	}
    return closestObject;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    Ray ray; // TODO: set up the initial view ray here
    ray.endpoint = camera.position;
    //vec3 w = camera.World_Position(pixel_index);
    ray.direction = (camera.World_Position(pixel_index) - camera.position).normalized();
    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    // TODO
    Hit temp;
    Object* foundobj = Closest_Intersection(ray, temp);
    vec3 color;
    
    if(foundobj != nullptr && recursion_depth <= recursion_depth_limit){
      vec3 intersection_point = ray.Point(temp.t);
      vec3 norm = foundobj->Normal(intersection_point);
      if(temp.ray_exiting){
        norm = norm * -1;
        }
      
      color = foundobj->material_shader->Shade_Surface(ray,intersection_point,norm,recursion_depth,temp.ray_exiting);
    }
    else{
      vec3 dummy;
      color = background_shader->Shade_Surface(ray,dummy,dummy,recursion_depth, temp.ray_exiting);
    }
    // determine the color here

    return color;
}
