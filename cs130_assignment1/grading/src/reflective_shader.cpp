#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth,bool is_exiting) const
{
    //r = v-2(dotvn)n
    vec3 color, v, reflected_color, shader_color, n;
    // TODO: determine the color
    
    v = ray.direction.normalized();
    n = same_side_normal.normalized();
    
    Ray r;
    r.endpoint = intersection_point;
    r.direction = v - 2.0 * dot(v, n) * n;
    //++recursion_depth;
    reflected_color = world.Cast_Ray(r, recursion_depth + 1);
    shader_color = shader->Shade_Surface(ray, intersection_point, same_side_normal, recursion_depth, is_exiting);
    color = reflectivity * reflected_color + (1.0 - reflectivity) * shader_color;
  

   
    
    return color;
}
