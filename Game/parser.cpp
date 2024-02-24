#include <fstream>
#include "parser.h"
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>

using namespace std;

class parser {
public:
    Eye* eye;
    vec4* ambientLight;
    vector<Light*>* lights;
    vector<SpotLight*>* spotlights;
    vector<Sphere*>* spheres;
    vector<Plane*>* planes;
    vector<Surface*>* objects = new vector<Surface*>();

     parser() {
        this->eye = new Eye();
        this->ambientLight = new glm::vec4(0.0f);
        this->lights = new vector<Light*>();
        this->spotlights = new vector<SpotLight*>();
        this->spheres = new vector<Sphere*>();
        this->planes = new vector<Plane*>();
    };
 
      void parseInputFile(string fileName)
      {
          int colindex = 0, posindex = 0, intindex = 0;
        

          char l = 'a';
          double x1 = 1;
          double x2 = 1;
          double x3 = 1;
          double x4 = 1;
          std::ifstream inputFile(fileName); // Open the file for reading
          if (inputFile.is_open()) {
              std::cout << "its open" << std::endl;
          }
          if (!inputFile) {
              std::cerr << "Error opening file " << fileName << ": " << std::strerror(errno) << std::endl;
              
          }

          std::string line;
          while (std::getline(inputFile, line)) {
              std::cout << line << std::endl;
              char firstChar = 'a';
              vector<double> numbers; // Store numbers in a vector

              std::istringstream iss(line);
              iss >> firstChar; // Extract the first character

              double number;
              while (iss >> number) {
                  numbers.push_back(number); // Store the number in the vector
              }

              // Output the extracted character and numbers
              std::cout << "First character: " << firstChar << std::endl;
              std::cout << "Numbers: ";
              std::vector<double>::iterator it;
              int i = 0;
              for (it = numbers.begin(); it != numbers.end(); ++it) {
                  std::cout << numbers[i] << " ";
                  i++;
              }

              std::cout << std::endl;
              x1 = numbers[0];
              x2 = numbers[1];
              x3 = numbers[2];
              x4 = numbers[3];

              /*std::cout << "the one is " << one << "/n" << std::endl;
              std::cout << "the two is " << two << "/n" << std::endl;
              std::cout << "the three is " << three << "/n" << std::endl;
              std::cout << "the four is " << four << "/n" << std::endl;*/
              l = firstChar;
          
              switch (l)
              {
              case 'e':
                  this->eye = parseEye(x1, x2, x3, x4);
                  break;
              case 'a':
                  this->ambientLight = new vec4(x1, x2, x3, x4);
                  break;
              case 'd':
                  parseLightDirection(x1, x2, x3, x4);

                  break;
              case 'p':
                  this->spotlights->at(posindex)->setPosition(x1, x2, x3);
                  this->spotlights->at(posindex)->setAngle(x4);
                  (posindex)++;
                  break;
              case 'i':
                  this->lights->at(intindex)->setIntensity(vec4(x1, x2, x3, x4));
                  (intindex)++;
                  break;
              case 'c':
                  this->objects->at(colindex)->setColor(vec4(x1, x2, x3, x4));
                  this->objects->at(colindex)->setShininess(x4);
                  (colindex)++;
                  break;
              default:
                  //objectDescriptor obj;
                  if (x4 > 0)
                  {
                      Sphere* obj = parseSphere(l, x1, x2, x3, x4);
                      this->spheres->push_back(obj);
                      obj = this->spheres->at(this->spheres->size() - 1);
                      this->objects->push_back(obj);
                      //  
                  }
                  else
                  {
                      Plane* obj = parsePlane(l, x1, x2, x3, x4);
                      this->planes->push_back(obj);

                      obj = this->planes->at(this->planes->size() - 1);
                      this->objects->push_back(obj);
                      // 
                  }
              }
          }
      }

      void parseLine(char op, double x1, double x2, double x3, double x4,
          vector<Surface*>* orderedObjects, int* nextColoredObject, int* nextSpotlightPosition, int* nextIntensityPosition)
      {

          
      }

      static Eye* parseEye(double x1, double x2, double x3, double x4)
      {
          return new Eye(x1, x2, x3);
      }

      Light* parseLightDirection(double x1, double x2, double x3, double x4)
      {
          Light* light = nullptr;
          if (x4 == 1) {
              light = new SpotLight(vec3(x1, x2, x3));
              this->spotlights->push_back(const_cast<SpotLight*>(reinterpret_cast<const SpotLight*>(light)));
          }
          else
              light = new DirectionalLight(vec3(x1,x2,x3));
          light->setDirection(x1, x2, x3);
          this->lights->push_back(light);
          return light;
      }

      static ObjectType getType(char c)
      {
          switch (c)
          {
          case 'r':
              return Reflective;
          case 't':
              return Transparent;
          case 'o':
              return Object;
          default:
              return Nothing;
          }
      }

      static Sphere* parseSphere(char op, double x1, double x2, double x3, double x4)
      {
          ObjectType type = getType(op);
          Sphere* s = new Sphere(x1, x2, x3, x4, type);
          s->setRadius(x4);
          return s;
      }

      static Plane* parsePlane(char op, double x1, double x2, double x3, double x4)
      {
          ObjectType type = getType(op);
          Plane* p = new Plane(x1, x2, x3, x4, type);
          return p;
      }

};