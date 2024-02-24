#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <glm/glm.hpp>
#include <string>
#pragma once

using namespace std;
using namespace glm;

enum ObjectType { Object, Reflective, Transparent , Nothing};
enum ObjectClass { PlaneObject, SphereObject };
enum LightType { Directional, Spotlight };

struct Surface
{
protected:
    ObjectType type;
    ObjectClass objectClass = PlaneObject;
    vec4 coordinates = vec4(0, 0, 0, 0);
    vec3 color = vec3(0, 0, 0);
    vec3 position = vec3(0, 0, 0);
    float shininess = 0;

public:
    
    void setColor(vec4 color) {
        this->color = vec3(color.r, color.g, color.b);
        this->shininess = color.w;
    }
    void setShininess(float shininess) {
        this->shininess = shininess;
    }
    float getShininess() {
        return this->shininess;
    }
    vec4 getCoordinates() {
        return this->coordinates;
    }
    ObjectType getType() {
        return this->type;
    }
    ObjectClass getClassOfObject() {
        return this->objectClass;
    }
    vec3 getPosition() {
        return position;
    }
    virtual vec3 getColor(vec3 hit_point) = 0;
};


struct Sphere : Surface
{
private:
    double radius;

public:
    Sphere(double x, double y, double z, double r, ObjectType type) {
        this->coordinates = vec4(x, y, z, r);
        this->position = glm::vec3(x, y, z);
        this->radius = r;
        this->type = type;
        this->objectClass = SphereObject;
    };

    void setRadius(double r){ 
        this->radius = r;
    }

    vec3 getPosition() {
        return this->position;
    }

    float getRadius() {
        return this->radius;
    }

    vec3 getColor(vec3 hitPoint) {
        return this->color;
    }

};

struct Plane : Surface
{

    Plane(double a, double b, double c, double d, ObjectType type) {
        this->coordinates = vec4(a, b, c, d);
        this->position = vec3(a, b, c);
        this->type = type;
        this->objectClass = PlaneObject;
    };

    vec3 getPosition() {
        return this->position;
    }
    float getD() {
        return this->coordinates.w;
    }
    ////////////from the tergol 5 /////////////////
    vec3 getColor(vec3 hitPoint) {
        // Checkerboard pattern
        float scale_parameter = 0.5f;
        float chessboard = 0;

        if (hitPoint.x < 0) {
            chessboard += floor((0.5 - hitPoint.x) / scale_parameter);
        }
        else {
            chessboard += floor(hitPoint.x / scale_parameter);
        }

        if (hitPoint.y < 0) {
            chessboard += floor((0.5 - hitPoint.y) / scale_parameter);
        }
        else {
            chessboard += floor(hitPoint.y / scale_parameter);
        }

        chessboard = (chessboard * 0.5) - int(chessboard * 0.5);
        chessboard *= 2;
        if (chessboard > 0.5) {
            return 0.5f * this->color;
        }
        return this->color;
    }
  
};



struct Ray {

private:
    vec3 direction;
    vec3 origin;
    vec3 hitPoint;
    Surface* sceneObject;

public:
    Ray(vec3 direction, vec3 origin) {
        this->direction = direction;
        this->origin = origin;
        this->hitPoint = origin + direction;
        this->sceneObject = new Plane(0.0, 0.0, 0.0, 0.0, Nothing);;
    }

    vec3 getRayDirection() {
        return this->direction;
    }
    vec3 getRayOrigin() {
        return this->origin;
    }
    vec3 getHitPoint() {
        return this->hitPoint;
    }
    Surface* getSceneObject() {
        return this->sceneObject;
    }
    void setRayDirection(vec3 vec) {
        this->direction = vec;
    }

    void setRayOrigin(vec3 vec) {
        this->origin = vec;
    }
    void setHitPoint(vec3 vec) {
        this->hitPoint = vec;
    }
    void setSceneObject(Surface* obj) {
        this->sceneObject = obj;
    }
    
};


struct Eye 
{
    vec3 coordinates;

    Eye(double x, double y, double z) {
        coordinates = vec3(x, y, z);
    };
    Eye() {
        coordinates = vec3(0, 0, 0);
    };
    vec3 getCoordinates() {
       return this->coordinates;
    }


};

struct Light 
{
    glm::vec3 direction;
    glm::vec3 intensity;
    float shininess = 0;
    LightType type;

    void setDirection(float x, float y, float z) {
        this->direction = glm::vec3(x, y, z);
    }

    void setIntensity(vec4 intensity) {
        this->intensity = vec3(intensity.r, intensity.g, intensity.b);
        this->shininess = intensity.w;
    }

    vec3 getIntensity() {
        return this->intensity;
    }
};

struct DirectionalLight : Light
{
    DirectionalLight(vec3 direction) {
        this->type = Directional;
        this->direction = direction;
    }
};

struct SpotLight : Light
{
    vec3 position;
    float w = 0;

    SpotLight(vec3 direction) {
        this->type = Spotlight;
        this->direction = direction;
        w = 0;
        position = vec3(0, 0, 0);
    }
    

    void setPosition(float x, float y, float z) {
        this->position = glm::vec3(x, y, z);
    }
    void setAngle(float w) {
        this->w = w;
    }
    float getAngle() {
        return w;
    }
    vec3 getPosition() {
        return this->position;
    }
};


#endif
