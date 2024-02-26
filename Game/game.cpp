#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "parser.h"
#include "parser.cpp"
#include "algorithm"

Ray UpdateRay(int j, int i, Surface* ob, bool update, Ray reflectedRay, parser* scene) {

    float width = 2.0f / 800.0f;
    float height = 2.0f / 800.0f;

    if (!update) {
        vec3 pixelCenter(-1 + width / 2, 1 - height / 2, 0);
        vec3 exactPixel = pixelCenter + vec3(j * width, -1 * (i * height), 0);
        vec3 eyeVec = scene->eye->getCoordinates();
        vec3 rayDirection = normalize(exactPixel - eyeVec);

        reflectedRay.setRayDirection(rayDirection);
        reflectedRay.setRayOrigin(eyeVec);

    }

    // update the ray
    float minT = INFINITY;
    Surface* closestObject = new Plane(0.0, 0.0, 0.0, 0.0, Nothing);
    reflectedRay.setHitPoint(reflectedRay.getRayOrigin() + reflectedRay.getRayDirection());
    reflectedRay.setSceneObject(closestObject);

    for (int i = 0; i < scene->objects->size(); i++) {
        float t = 0.0;
        Surface* currentObject = scene->objects->at(i);
        if (currentObject != ob) {

            if (currentObject->getClassOfObject() == PlaneObject) {
                float denominator = glm::dot(reflectedRay.getRayDirection(), currentObject->getPosition());

                if (abs(denominator) < 0.0001f) {
                    t = -1.0f; // No intersection

                }
                // intersection equation
                t = -(glm::dot(reflectedRay.getRayOrigin(), currentObject->getPosition()) + ((Plane*)currentObject)->getD()) / denominator;

                if (t < 0.0f) {
                    t = -1.0f; // No intersection

                }

            }

            else {
                vec3 oc = reflectedRay.getRayOrigin() - currentObject->getPosition();
                float a = dot(reflectedRay.getRayDirection(), reflectedRay.getRayDirection());
                float b = 2.0f * dot(oc, reflectedRay.getRayDirection());
                float c = dot(oc, oc) - ((Sphere*)currentObject)->getRadius() * ((Sphere*)currentObject)->getRadius();

                float delta = b * b - 4 * a * c; // Discriminant of the quadratic equation

                if (delta < 0) {
                    t = -1.0f;

                }
                else {
                    float t1 = (-b - sqrt(delta)) / (2.0f * a);
                    float t2 = (-b + sqrt(delta)) / (2.0f * a);

                    if (t1 < 0 && t2 < 0) {
                        t = -1.0f; // no intersection

                    }
                    float first = (t1 >= 0) ? t1 : t2;

                    if (first <= 0.0001f) {
                        float second = (t1 >= 0 && t2 >= 0) ? glm::max(t1, t2) : -1.0f;
                        t = second;

                    }
                    else {
                        t = first;

                    }
                }
            }
            if ((t >= 0) && t < minT) {
                closestObject = currentObject;
                minT = t;
                reflectedRay.setSceneObject(closestObject);
                reflectedRay.setHitPoint(reflectedRay.getRayOrigin() + reflectedRay.getRayDirection() * minT);
            }
        }
    }

    return reflectedRay;
}

vec3 getNormal(vec3 hit_point, Surface* obj) {
    if (obj->getClassOfObject() == PlaneObject) {
        return normalize(vec3(obj->getCoordinates()));
    }
    else
        return normalize(hit_point - ((Sphere*)obj)->getPosition());
}

float calcNLi(vec3 N, Ray ray, Light* light) { //defuse

    if (ray.getSceneObject()->getClassOfObject() == PlaneObject) {
        vec3 lightDirection = -normalize(light->direction);
        if (light->type == Directional) {
            float cos = dot(N, -lightDirection);
            return glm::max(cos, 0.0f);
        }
        else {
            vec3 SpotRay = normalize(ray.getHitPoint() - ((SpotLight*)light)->getPosition());
            float cos = dot(SpotRay, -lightDirection);
            if (cos < ((SpotLight*)light)->getAngle())
                return 0.0f;
            else {
                lightDirection = -SpotRay;
                cos = dot(N, -lightDirection);
                return glm::max(cos, 0.0f);;
            }

        }
    }
    else { // it is a sphere
        vec3 lightDirection = normalize(light->direction);
        if (light->type == Directional) {
            float cos = dot(N, -lightDirection);
            return glm::max(cos, 0.0f);;
        }
        else {
            vec3 SpotRay = normalize(ray.getHitPoint() - ((SpotLight*)light)->getPosition());
            float cos = dot(SpotRay, lightDirection);
            if (cos < ((SpotLight*)light)->getAngle())
                return 0.0f;
            else {
                lightDirection = SpotRay;
                cos = dot(N, -lightDirection);
                return glm::max(cos, 0.0f);;
            }
        }
    }
}

float calcVRi(vec3 V, Ray ray, Light* light) { //specular

    if (ray.getSceneObject()->getClassOfObject() == PlaneObject) { // it is a plane
        vec3 lightDirection = normalize(light->direction);
        vec3 normalPlane = getNormal(ray.getHitPoint(), ray.getSceneObject());
        if (light->type == Directional) {
            vec3 reflectedEquation = lightDirection - 2.0f * normalPlane * dot(lightDirection, normalPlane);
            float cos = dot(V, reflectedEquation);
            cos = glm::max(0.0f, cos);
            cos = pow(cos, ray.getSceneObject()->getShininess());
            return cos;
        }
        else {
            vec3 SpotRay = normalize(ray.getHitPoint() - ((SpotLight*)light)->getPosition());
            float cos = dot(SpotRay, lightDirection);
            if (cos < ((SpotLight*)light)->getAngle()) {
                return 0.0f;
            }
            else {
                lightDirection = SpotRay;
                vec3 reflectedEquation = lightDirection - 2.0f * normalPlane * dot(lightDirection, normalPlane);
                float cos = dot(V, reflectedEquation);
                cos = glm::max(0.0f, cos);
                cos = pow(cos, ray.getSceneObject()->getShininess());
                return cos;
            }
        }
    }
    else { // it is a sphere
        vec3 lightDirection = normalize(light->direction);
        vec3 normalSphere = getNormal(ray.getHitPoint(), ray.getSceneObject());
        if (light->type == Directional) {
            vec3 reflectedEquation = lightDirection - 2.0f * normalSphere * dot(lightDirection, normalSphere);
            float cos = dot(V, reflectedEquation);
            cos = glm::max(0.0f, cos);
            cos = pow(cos, ray.getSceneObject()->getShininess());
            return cos;
        }
        else {
            vec3 SpotRay = normalize(ray.getHitPoint() - ((SpotLight*)light)->getPosition());
            float cos = dot(SpotRay, lightDirection);
            if (cos < ((SpotLight*)light)->getAngle()) {
                return 0.0f;
            }
            else {
                lightDirection = SpotRay;
                vec3 reflectedEquation = lightDirection - 2.0f * normalSphere * dot(lightDirection, normalSphere);
                float cos = dot(V, reflectedEquation);
                cos = glm::max(0.0f, cos);
                cos = pow(cos, ray.getSceneObject()->getShininess());
                return cos;
            }
        }
    }
}

float calcSiIi(Ray ray, Light* light, parser* scene) { //shadow

    vec3 lightDirection = glm::normalize(light->direction);
    float minT = INFINITY;

    if (light->type == Spotlight) {
        vec3 SpotRay = glm::normalize(ray.getHitPoint() - ((SpotLight*)light)->getPosition());
        float cos = dot(SpotRay, lightDirection);

        if (cos < ((SpotLight*)light)->getAngle()) {
            return 0.0;
        }
        else {
            lightDirection = SpotRay;
            minT = glm::length(((SpotLight*)light)->getPosition() - ray.getHitPoint());
        }
    }

    for (int i = 0; i < scene->objects->size(); i++) {
        Surface* currentObject = scene->objects->at(i);

        if (currentObject != ray.getSceneObject()) {
            Ray ray2 = Ray(-lightDirection, ray.getHitPoint());
            //float t = Intersection(scene->objects->at(i), ray2);
            float t = 0.0;

            if (currentObject != ray.getSceneObject()) {

                if (currentObject->getClassOfObject() == PlaneObject) {
                    float denominator = glm::dot(ray2.getRayDirection(), currentObject->getPosition());

                    if (abs(denominator) < 0.0001f) {
                        t = -1.0f; // No intersection

                    }
                    // intersection equation
                    t = -(glm::dot(ray2.getRayOrigin(), currentObject->getPosition()) + ((Plane*)currentObject)->getD()) / denominator;

                    if (t < 0.0f) {
                        t = -1.0f; // No intersection

                    }

                }

                else {
                    vec3 oc = ray2.getRayOrigin() - currentObject->getPosition();
                    float a = dot(ray2.getRayDirection(), ray2.getRayDirection());
                    float b = 2.0f * dot(oc, ray2.getRayDirection());
                    float c = dot(oc, oc) - ((Sphere*)currentObject)->getRadius() * ((Sphere*)currentObject)->getRadius();

                    float delta = b * b - 4 * a * c; // Discriminant of the quadratic equation

                    if (delta < 0) {
                        t = -1.0f;

                    }
                    else {
                        float t1 = (-b - sqrt(delta)) / (2.0f * a);
                        float t2 = (-b + sqrt(delta)) / (2.0f * a);

                        if (t1 < 0 && t2 < 0) {
                            t = -1.0f; // no intersection

                        }
                        float first = (t1 >= 0) ? t1 : t2;

                        if (first <= 0.0001f) {
                            float second = (t1 >= 0 && t2 >= 0) ? glm::max(t1, t2) : -1.0f;
                            t = second;

                        }
                        else {
                            t = first;

                        }
                    }
                }
            }

            if ((t > 0) && (t < minT)) {
                return 0.0;
            }

        }
    }

    return 1.0;
}

Ray SnellLaw(Ray ray, glm::vec3 N, glm::vec3 rayDirection, float snellFrac) {
    const float PI = 3.14159265;

    // Snell's Law calculations
    vec3 surfaceNormal = getNormal(ray.getHitPoint(), ray.getSceneObject());
    float cos1 = dot(surfaceNormal, -ray.getRayDirection());
    float theta1 = acos(cos1) * (180.0f / PI);
    float snellFraction = (1.0f / 1.5f);
    float sin1 = sin(theta1 * (PI / 180.0f));
    float sin2 = snellFraction * sin1;
    float theta2 = asin(sin2) * (180.0f / PI);
    float cos2 = cos(theta2 * (PI / 180.0f));

    vec3 direction = (snellFrac * cos1 - cos2) * N - snellFrac * (-ray.getRayDirection());
    vec3 org = ray.getHitPoint();
    Ray newRay(direction, org);
    return newRay;
}

vec4 GetPixelColor(int j, int i, Ray ray, int count, parser* scene) {
    vec3 I(0, 0, 0);
    vec3 IE(0, 0, 0);
    vec3 KA(0, 0, 0);
    vec3 IA(0, 0, 0);
    vec3 KR(0, 0, 0);
    vec3 IR(0, 0, 0);
    vec3 X1(0, 0, 0); //N*L
    vec3 X2(0, 0, 0); //(V*R)^n
    vec3 segma(0, 0, 0);
    if (ray.getSceneObject()->getType() == Object) { // Object Type
        KA = ray.getSceneObject()->getColor(ray.getHitPoint());
        IA = vec3(scene->ambientLight->r, scene->ambientLight->g, scene->ambientLight->b);

        for (int i = 0; i < scene->lights->size(); i++) {
            vec3 KS(0.7f, 0.7f, 0.7f);
            vec3 KD = ray.getSceneObject()->getColor(ray.getHitPoint()) * scene->lights->at(i)->getIntensity();
            KS = KS * scene->lights->at(i)->getIntensity();
            vec3 N = getNormal(ray.getHitPoint(), ray.getSceneObject());
            vec3 V = normalize(ray.getRayOrigin() - ray.getHitPoint());
            X1 = KD * calcNLi(N, ray, scene->lights->at(i));
            X2 = (KS * calcVRi(V, ray, scene->lights->at(i)));
            float X3 = calcSiIi(ray, scene->lights->at(i), scene);
            segma = segma + (X1 + X2) * X3;
        }
    }
    I = IE + (KA * IA) + segma + (KR * IR);

    if (ray.getSceneObject()->getType() == Reflective) { // Reflective Type
        if (count == 5) {
            return vec4(0.f, 0.f, 0.f, 0.f);
        }
        // update the ray
        vec3 reflectionEquation = ray.getRayDirection() - 2.0f * getNormal(ray.getHitPoint(), ray.getSceneObject()) * dot(ray.getRayDirection(), getNormal(ray.getHitPoint(), ray.getSceneObject()));
        Ray newRay(reflectionEquation, ray.getHitPoint());
        newRay = UpdateRay(j, i, ray.getSceneObject(), true, newRay, scene);

        if (newRay.getSceneObject()->getType() == Nothing) {
            return vec4(0.f, 0.f, 0.f, 0.f);
        }

        vec4 newColor = GetPixelColor(j, i, newRay, count + 1, scene);
        I = vec3(newColor.r, newColor.g, newColor.b);
    }

    if (ray.getSceneObject()->getType() == Transparent) {// Transparent Type

        vec3 surfaceNormal = getNormal(ray.getHitPoint(), ray.getSceneObject());
        float snellFraction = (1.0f / 1.5f);
        Ray insideRay = SnellLaw(ray, surfaceNormal, ray.getRayDirection(), snellFraction);

        Plane* blackPlane = new Plane(0.0, 0.0, 0.0, 0.0, Nothing);
        insideRay = UpdateRay(j, i, blackPlane, true, insideRay, scene);
        Surface* obj = ray.getSceneObject();
        float t = 0.0;
        if (obj->getClassOfObject() == PlaneObject) {
            float denominator = glm::dot(ray.getRayDirection(), obj->getPosition());

            if (abs(denominator) < 0.0001f) {
                t = -1.0f; // No intersection
            }
            float t = -(glm::dot(ray.getRayOrigin(), obj->getPosition()) + ((Plane*)obj)->getD()) / denominator;

            if (t < 0.0f) {
                t = -1.0f; // No intersection
            }
        }

        else {
            vec3 oc = ray.getRayOrigin() - obj->getPosition();
            float a = dot(ray.getRayDirection(), ray.getRayDirection());
            float b = 2.0f * dot(oc, ray.getRayDirection());
            float c = dot(oc, oc) - ((Sphere*)obj)->getRadius() * ((Sphere*)obj)->getRadius();

            float delta = b * b - 4 * a * c; // Discriminant of the quadratic equation

            if (delta < 0) {
                t = -1.0f;
            }
            else {
                float t1 = (-b - sqrt(delta)) / (2.0f * a);
                float t2 = (-b + sqrt(delta)) / (2.0f * a);

                if (t1 < 0 && t2 < 0) {
                    t = -1.0f; // no intersection
                }
                float first = (t1 >= 0) ? t1 : t2;

                if (first <= 0.0001f) {
                    float second = (t1 >= 0 && t2 >= 0) ? glm::max(t1, t2) : -1.0f;
                    t = second;
                }
                else {
                    t = first;
                }
            }
        }
        vec3 secondHit = insideRay.getRayOrigin() + insideRay.getRayDirection() * t;
        surfaceNormal = getNormal(secondHit, insideRay.getSceneObject());
        Ray refractedRay = SnellLaw(insideRay, -surfaceNormal, insideRay.getRayDirection(), snellFraction);

        refractedRay = UpdateRay(j, i, insideRay.getSceneObject(), true, insideRay, scene);
        if (refractedRay.getSceneObject()->getType() == Nothing) {
            return vec4(0.f, 0.f, 0.f, 0.f);
        }
        vec4 color = GetPixelColor(j, i, refractedRay, count + 1, scene);


        I = vec3(color.r, color.g, color.b);
    }

    I = min(I, vec3(1.0, 1.0, 1.0));
    I = max(I, vec3(0.0, 0.0, 0.0));
    return vec4(I.r, I.g, I.b, 1.0);


}
unsigned char* rendering(parser* scene) {
    auto* Data = new unsigned char[800 * 800 * 4];

    for (int i = 0; i < 800; i++) {
        for (int j = 0; j < 800; j++) {

            Ray initial_ray(vec3(0, 0, 0), vec3(0, 0, 0));
            Plane* blackPlane = new Plane(0.0, 0.0, 0.0, 0.0, Nothing);
            Ray ray = UpdateRay(j, i, blackPlane, false, initial_ray, scene);
            vec4 color = GetPixelColor(j, i, ray, 0, scene);

            Data[(j + 800 * i) * 4] = (unsigned char)(color.r * 255);
            Data[(j + 800 * i) * 4 + 1] = (unsigned char)(color.g * 255);
            Data[(j + 800 * i) * 4 + 2] = (unsigned char)(color.b * 255);
            Data[(j + 800 * i) * 4 + 3] = (unsigned char)(color.a * 255);

        }

    }
    return Data;
}

static void printMat(const glm::mat4 mat)
{
    std::cout << " matrix:" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            std::cout << mat[j][i] << " ";
        std::cout << std::endl;
    }
}

Game::Game() : Scene()
{
}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1)
{
}

void Game::Init()
{

    AddShader("../res/shaders/pickingShader");
    AddShader("../res/shaders/basicShader");
    parser* p = new parser();
    p->parse("../res/Scenes/scene4.txt");
    unsigned char* data = rendering(p);

    AddTexture(800, 800, data);
    AddShape(Plane, -1, TRIANGLES);

    pickedShape = 0;

    SetShapeTex(0, 0);
    MoveCamera(0, zTranslate, 10);
    pickedShape = -1;

    //ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4& MVP, const glm::mat4& Model, const int  shaderIndx)
{
    Shader* s = shaders[shaderIndx];
    int r = ((pickedShape + 1) & 0x000000FF) >> 0;
    int g = ((pickedShape + 1) & 0x0000FF00) >> 8;
    int b = ((pickedShape + 1) & 0x00FF0000) >> 16;
    s->Bind();
    s->SetUniformMat4f("MVP", MVP);
    s->SetUniformMat4f("Normal", Model);
    s->SetUniform4f("lightDirection", 0.0f, 0.0f, -1.0f, 0.0f);
    if (shaderIndx == 0)
        s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
    else
        s->SetUniform4f("lightColor", 0.7f, 0.8f, 0.1f, 1.0f);
    s->Unbind();
}

void Game::WhenRotate()
{
}

void Game::WhenTranslate()
{
}

void Game::Motion()
{
    if (isActive)
    {
    }
}

Game::~Game(void)
{
}
