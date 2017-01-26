/*---------------------------------------------------------------------
*
* Copyright © 2015  Minsi Chen
* E-mail: m.chen@derby.ac.uk
*
* The source is written for the Graphics I and II modules. You are free
* to use and extend the functionality. The code provided here is functional
* however the author does not guarantee its performance.
---------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#if defined(WIN32) || defined(_WINDOWS)
#include <Windows.h>
#include <gl/GL.h>
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#endif

#include "RayTracer.h"
#include "Ray.h"
#include "Scene.h"
#include "Camera.h"
#include "perlin.h"
#define PI 3.14159265

RayTracer::RayTracer()
{
	m_buffHeight = m_buffWidth = 0.0;
	m_renderCount = 0;
	SetTraceLevel(5);
	m_traceflag = (TraceFlags)(TRACE_AMBIENT | TRACE_DIFFUSE_AND_SPEC |
		TRACE_SHADOW | TRACE_REFLECTION | TRACE_REFRACTION);

}

RayTracer::RayTracer(int Width, int Height)
{
	m_buffWidth = Width;
	m_buffHeight = Height;
	m_renderCount = 0;
	SetTraceLevel(5);

	m_framebuffer = new Framebuffer(Width, Height);

	//default set default trace flag, i.e. no lighting, non-recursive
	m_traceflag = (TraceFlags)(TRACE_AMBIENT);
}

RayTracer::~RayTracer()
{
	delete m_framebuffer;
}

void RayTracer::DoRayTrace(Scene* pScene)
{
	Camera* cam = pScene->GetSceneCamera();

	Vector3 camRightVector = cam->GetRightVector();
	Vector3 camUpVector = cam->GetUpVector();
	Vector3 camViewVector = cam->GetViewVector();
	Vector3 centre = cam->GetViewCentre();
	Vector3 camPosition = cam->GetPosition();

	double sceneWidth = pScene->GetSceneWidth();
	double sceneHeight = pScene->GetSceneHeight();

	double pixelDX = sceneWidth / m_buffWidth;
	double pixelDY = sceneHeight / m_buffHeight;

	int total = m_buffHeight*m_buffWidth;
	int done_count = 0;

	Vector3 start;

	start[0] = centre[0] - ((sceneWidth * camRightVector[0])
		+ (sceneHeight * camUpVector[0])) / 2.0;
	start[1] = centre[1] - ((sceneWidth * camRightVector[1])
		+ (sceneHeight * camUpVector[1])) / 2.0;
	start[2] = centre[2] - ((sceneWidth * camRightVector[2])
		+ (sceneHeight * camUpVector[2])) / 2.0;

	Colour scenebg = pScene->GetBackgroundColour();

	if (m_renderCount == 0)
	{
		fprintf(stdout, "Trace start.\n");

		Colour colour;
		//TinyRay on multiprocessors using OpenMP!!!
#pragma omp parallel for schedule (dynamic, 1) private(colour)
		for (int i = 0; i < m_buffHeight; i += 1) {
			for (int j = 0; j < m_buffWidth; j += 1) {

				//calculate the metric size of a pixel in the view plane (e.g. framebuffer)
				Vector3 pixel;

				pixel[0] = start[0] + (i + 0.5) * camUpVector[0] * pixelDY
					+ (j + 0.5) * camRightVector[0] * pixelDX;
				pixel[1] = start[1] + (i + 0.5) * camUpVector[1] * pixelDY
					+ (j + 0.5) * camRightVector[1] * pixelDX;
				pixel[2] = start[2] + (i + 0.5) * camUpVector[2] * pixelDY
					+ (j + 0.5) * camRightVector[2] * pixelDX;

				/*
				* setup first generation view ray
				* In perspective projection, each view ray originates from the eye (camera) position
				* and pierces through a pixel in the view plane
				*/
				Ray viewray;
				viewray.SetRay(camPosition, (pixel - camPosition).Normalise());

				double u = (double)j / (double)m_buffWidth;
				double v = (double)i / (double)m_buffHeight;

				scenebg = pScene->GetBackgroundColour();

				//trace the scene using the view ray
				//default colour is the background colour, unless something is hit along the way
				colour = this->TraceScene(pScene, viewray, scenebg, m_traceLevel);

				/*
				* Draw the pixel as a coloured rectangle
				*/
				m_framebuffer->WriteRGBToFramebuffer(colour, j, i);
			}
		}

		fprintf(stdout, "Done!!!\n");
		m_renderCount++;
	}
}

Colour RayTracer::TraceScene(Scene* pScene, Ray& ray, Colour incolour, int tracelevel, bool shadowray)
{
	if (tracelevel <= 0) // having 5 trace level we descrease them everytime we call the function if it hits 0 it returns  
		return incolour;

	RayHitResult result;

	Colour outcolour = incolour; //the output colour based on the ray-primitive intersection

	std::vector<Light*> *light_list = pScene->GetLightList();
	std::vector<Light*>::iterator lit_iter = light_list->begin();
	//Vector3 cameraPosition = pScene->GetSceneCamera()->GetPosition(); 

	//Intersect the ray with the scene
	result = pScene->IntersectByRay(ray);

	if (result.data) //the ray has hit something
	{
		if (shadowray) //if it is a shadowray, the color is darken and it returns from the function
		{
			//if my shadowray hit something I need to check if the object hit can cast shadows
			if (((Primitive*)result.data)->GetMaterial()->CastShadow())
			{
				outcolour = outcolour * 0.3f;
			}
			return outcolour ;
		}

		//	 When a ray intersect with an objects, determine the colour at the intersection point
		//	 using CalculateLighting
		outcolour = CalculateLighting(light_list, &pScene->GetSceneCamera()->GetPosition(), &result);

		//Reflection
		if (m_traceflag & TRACE_REFLECTION)
		{
			if (((Primitive*)result.data)->m_primtype != Primitive::PRIMTYPE_Plane)
			{
				Vector3 reflection = ray.GetRay().Reflect(result.normal); // calculating reflection vector
				Vector3 noiseReflection = reflection * 0.001;			  //making sure that when a ray is shot they dont overlap
				ray.SetRay(result.point + noiseReflection, reflection);
				outcolour = outcolour + TraceScene(pScene, ray, outcolour, tracelevel - 1, shadowray) * 0.3;  //calline trace line passing in the outcolor and decreasing the trace level 
			}
		}

		//Refraction
		if (m_traceflag & TRACE_REFRACTION)
		{
			//TODO: trace the refraction ray from the intersection point
			if (((Primitive*)result.data)->m_primtype != Primitive::PRIMTYPE_Plane)
			{
				float n1 = 1.0002772;	//air index
				float n2 = 1.523f;		//glass index
				float cos = ray.GetRay().DotProduct(result.normal);	
				//Snell law -> n1 * sin(angleI) = n2 * sin(angleR) , we can divide both side by n2 and obtain (n1 * sin(angleI)/ n2 = sin(angleR))
				float angleRefract = n1 * sin(cos) / n2;				
				if (angleRefract > 1)
					angleRefract = 1;
				if (angleRefract < 0)
					angleRefract = 0;
				Vector3 refraction = ray.GetRay().Refract(result.normal, sin(angleRefract));    // calculating refraction vector 
				Vector3 noiseRefraction = refraction * 0.001;									//making sure that when a ray is shot they dont overlap
				ray.SetRay(result.point + noiseRefraction, refraction);
				outcolour = outcolour + TraceScene(pScene, ray, outcolour, tracelevel - 1, shadowray) * 0.5;  //calline trace line passing in the outcolor and decreasing the trace level
			}
		}

		//Shadows
		if (m_traceflag & TRACE_SHADOW)
		{
			for each (Light* light in (*light_list)) 
			{
				//Calculating direction + normalisation
				Vector3 dir = (light->GetLightPosition() - result.point).Normalise(); 
				Vector3 shadowNoise = dir * 0.001;
				//Shooting a ray from the hit point to the light source to see if we intersect anything
				ray.SetRay(result.point + shadowNoise, dir);				
				//call TraceScene this time setting shadowray to true
				outcolour = TraceScene(pScene, ray, outcolour, tracelevel - 1, true);   
			}			
		}
	}

	return outcolour;
}

Colour RayTracer::CalculateLighting(std::vector<Light*>* lights, Vector3* campos, RayHitResult* hitresult)
{
	Colour outcolour;
	std::vector<Light*>::iterator lit_iter = lights->begin();

	Primitive* prim = (Primitive*)hitresult->data;
	Material* mat = prim->GetMaterial();

	//setting all the colors
	Colour diffColor = mat->GetDiffuseColour();
	Colour specColor = mat->GetSpecularColour();

	//setting vectors
	Vector3 lightVector = ((*lit_iter)->GetLightPosition() - hitresult->point).Normalise();			//lightVector 
	Vector3 viewVector = (*campos - hitresult->point).Normalise();									//viewVector

	outcolour = mat->GetAmbientColour();

	////Go through all lights in the scene
	for each (Light* light in (*lights))
	{
		if (m_traceflag & TRACE_DIFFUSE_AND_SPEC)
		{
			//diffuse using Lambertian model, for specular, I'm using Phong model
			float angleDiff = lightVector.DotProduct(hitresult->normal);						    //Calculating the angle between the light vector and the normal of the hit point
			float angleSpec = viewVector.DotProduct((lightVector * -1).Reflect(hitresult->normal));	//Calculating the angle between the view vector and the reflected vector

			//clamp results 0 to 1 for both values
			if (angleDiff < 0)
				angleDiff = 0;
			if (angleDiff > 1)
				angleDiff = 1;
			if (angleSpec < 0)
				angleSpec = 0;
			if (angleSpec > 1)
				angleSpec = 1;

			float specValue = pow(angleSpec, mat->GetSpecPower());				  //Calculating specular power by doing cos	ⁿ where n is the spec power of the material 

			//Calculate all the colors 
			diffColor = diffColor * (light)->GetLightColour() * angleDiff;					//K * l * cos    // Diffuse Color
			specColor = specColor * (light)->GetLightColour() * angleSpec * specValue;      //K * l * cos	ⁿ  // Specular Color
			outcolour = outcolour + specColor + diffColor;									//adding colors together to get the final color
		}
	}

	//Generate the grid pattern on the plane
	if (((Primitive*)hitresult->data)->m_primtype == Primitive::PRIMTYPE_Plane)
	{
		int dx = hitresult->point[0] / 2.0;
		int dy = hitresult->point[1] / 2.0;
		int dz = hitresult->point[2] / 2.0;

		if (dx % 2 || dy % 2 || dz % 2)
		{
			outcolour = Vector3(0.1, 0.1, 0.1);
		}
		else
		{
			outcolour = mat->GetDiffuseColour();
		}
	}

	return outcolour;
}

