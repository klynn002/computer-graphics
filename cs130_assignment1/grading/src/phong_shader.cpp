#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth,bool is_exiting) const 
{
   vec3 color;
	color = world.ambient_color * world.ambient_intensity * this->color_ambient;   

	for(unsigned int i = 0; i < world.lights.size(); i++){
		vec3 l, n, lightColor, lightVector;	
		lightVector = (world.lights[i]->position - intersection_point).normalized();

		Hit hitShadow;
		Ray rayShadow(intersection_point, lightVector);

		rayShadow.endpoint = rayShadow.Point(.001);

		Object* obj = world.Closest_Intersection(rayShadow, hitShadow);

		if(!world.enable_shadows 
			|| (obj == NULL && world.enable_shadows) 	|| (((lightVector - intersection_point).magnitude() 
				<= (rayShadow.Point(hitShadow.t) - intersection_point).magnitude()) && world.enable_shadows)) {
			
			n = same_side_normal;
			l = world.lights.at(i)->position - intersection_point;
			lightColor = world.lights.at(i)->Emitted_Light(ray) / l.magnitude_squared();

			double dotln = dot(n, l.normalized());

			color += this->color_diffuse * lightColor * std::max(0.0, dotln);

		
			vec3 r = (2 * dot(n, l.normalized()) * n - l.normalized()).normalized();
      vec3 c = world.camera.position - intersection_point;
      double dotrc = dot(r, c.normalized());
      double spec_intens = pow(std::max(0.0, dotrc), specular_power);
			color += this->color_specular * lightColor * spec_intens;
	}
}


    
    return color;
}
