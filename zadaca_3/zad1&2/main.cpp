#include <iostream>
#include <limits>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include "algorithm"
#include "geometry.h"
#include "OBJParser.h"

using namespace std;

#define pi 3.14159265359

typedef std::vector<Vec3f> Image;

bool SMOOTH_SHADING = true;

float deg2rad(const float &deg)
{
    return deg * pi / 180;
}

struct Camera
{
    Vec3f pos;
    Vec3f dir;
    float rotation_angle;

    Vec3f u; // right
    Vec3f v; // up
    Vec3f w; // dir normalizirani
    Vec3f t; // normala na dir

    Camera(const Vec3f &pos, const Vec3f &dir, const float &rotation_angle) : pos(pos), dir(dir), rotation_angle(deg2rad(rotation_angle))
    {
        w = this->dir.normalize();

        t = Vec3f(1, 0, 0);

        u = cross(w, t);
        u.normalize();

        v = cross(u, w);
        v.normalize();
    }
};

struct Viewport
{
    unsigned int nx;
    unsigned int ny;
    float fov;
    Viewport(const unsigned int &nx, const unsigned int &ny, const float &fov) : nx(nx), ny(ny), fov(fov) {}
};

//https://freecutout.com/environment-maps/
//https://blog.michelanders.nl/search/label/ray%20tracing%20concepts
struct EnvironmentMap
{
    vector<vector<Vec3f>> image;
    float r;
    int Width, Height;

    EnvironmentMap(const char *filename, const float &r) : r(r), Width(Width), Height(Height)
    {
        string standard, width, height, maxcolor;
        ifstream file(filename, ifstream::binary);
        file >> standard >> width >> height >> maxcolor;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        Width = stoi(width);
        Height = stoi(height);

        vector<unsigned char> buffer(istreambuf_iterator<char>(file), {});

        for (int i = 0; i < Height; i++)
        {
            vector<Vec3f> row;
            for (int j = 0; j <= Width; j++)
            {
                float r = (float)buffer[i * (Width * 3) + (j * 3)] / 255;
                float g = (float)buffer[i * (Width * 3) + (j * 3) + 1] / 255;
                float b = (float)buffer[i * (Width * 3) + (j * 3) + 2] / 255;

                row.push_back(Vec3f(r, g, b));
            }
            image.push_back(row);
        }
        file.close();
    }

    Vec3f calc(const Vec3f &orig, const Vec3f &dir)
    {
        float u = 0.5 + 0.5 * atan2(dir.x, dir.z) / pi;
        float v = 0.5 - asin(dir.y) / pi;

        int i = floor(u * Width);
        int j = floor(v * Height);

        return image[j][i];
    }
};

struct Light
{
    Vec3f position;
    float intensity;

    Light(const Vec3f &position, const float &intensity) : position(position), intensity(intensity) {}
};
typedef std::vector<Light> Lights;

struct Material
{
    Vec2f albedo; // difuzni i spekularni koeficijenti refleksije
    Vec3f diffuse_color;
    float specular_exponent;
    float refraction_index;
    float alpha;

    Material(const Vec2f &a, const Vec3f &color, const float &coef, const float &refraction_index, const float &alpha) : albedo(a), diffuse_color(color), specular_exponent(coef), refraction_index(refraction_index), alpha(alpha) {}
    Material() : albedo(Vec2f(1, 0)), diffuse_color(), specular_exponent(1.f) {}
};

struct Object
{
    Material material;
    virtual bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t) = 0;
    virtual Vec3f normal(const Vec3f &p) const = 0;
};
typedef std::vector<Object *> Objects;

struct Sphere : Object
{
    Vec3f c; // centar
    float r; // radius

    Sphere(const Vec3f &c, const float &r, const Material &m) : c(c), r(r)
    {
        Object::material = m;
    }

    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t)
    {
        Vec3f v = c - p; // vektor izmedju izvora zrake i centra

        if (d * v < 0) // skalarni produkt
        {
            // sfera se nalazi iza zrake
            return false;
        }
        else
        {
            // izracunaj projekciju
            Vec3f pc = p + d * ((d * v) / d.norm());
            if ((c - pc) * (c - pc) > r * r)
            {
                // nema sjeciste
                return false;
            }
            else
            {
                float dist = sqrt(r * r - (c - pc) * (c - pc));

                if (v * v > r * r) // izvor pogleda izvan sfere
                {
                    t = (pc - p).norm() - dist;
                }
                else // izvor pogleda unutar sfere
                {
                    t = (pc - p).norm() + dist;
                }

                return true;
            }
        }
    }

    Vec3f normal(const Vec3f &p) const
    {
        return (p - c).normalize();
    }
};

struct Cuboid : Object
{
    Vec3f c1;
    Vec3f c2;

    Vec3f normal1;
    Vec3f normal2;
    Vec3f normal3;

    Cuboid(const Vec3f &c1, const Vec3f &c2, const Material &m) : c1(c1), c2(c2)
    {
        Object::material = m;
    }

    // p - pocetak zrake
    // d - direction
    // t - object distance
    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t)
    {
        float tNear = numeric_limits<float>::min();
        float tFar = numeric_limits<float>::max();

        float minX = min(c1.x, c2.x);
        float minY = min(c1.y, c2.y);
        float minZ = min(c1.z, c2.z);
        float maxX = max(c1.x, c2.x);
        float maxY = max(c1.y, c2.y);
        float maxZ = max(c1.z, c2.z);

        if (d.x == 0)
        {
            if (p.x < minX || p.x > maxX)
            {
                return false;
            }
        }
        else
        {
            float t1 = (minX - p.x) / d.x;
            float t2 = (maxX - p.x) / d.x;

            if (t1 > t2)
            {
                swap(t1, t2);
            }

            tNear = max(tNear, t1);
            tFar = min(tFar, t2);
            if (tNear > tFar || tFar < 0)
                return false;
        }

        t = tNear;

        if (d.y == 0)
        {
            if (p.y < minY || p.y > maxY)
            {
                return false;
            }
        }
        else
        {
            float t1 = (minY - p.y) / d.y;
            float t2 = (maxY - p.y) / d.y;

            if (t1 > t2)
            {
                swap(t1, t2);
            }

            tNear = max(tNear, t1);
            tFar = min(tFar, t2);
            if (tNear > tFar || tFar < 0)
                return false;
        }

        t = tNear;

        if (d.z == 0)
        {
            if (p.z < minZ || p.z > maxZ)
            {
                return false;
            }
        }
        else
        {
            float t1 = (minZ - p.z) / d.z;
            float t2 = (maxZ - p.z) / d.z;

            if (t1 > t2)
            {
                swap(t1, t2);
            }

            tNear = max(tNear, t1);
            tFar = min(tFar, t2);

            if (tNear > tFar || tFar < 0)
                return false;
        }

        t = tNear;

        return true;
    }

    Vec3f normal(const Vec3f &p) const
    {
        Vec3f normal;

        if (abs(p.x - c1.x) < 0.01)
            normal = Vec3f(-1, 0, 0);
        else if (abs(p.x - c2.x) < 0.01)
            normal = Vec3f(1, 0, 0);
        else if (abs(p.y - c1.y) < 0.01)
            normal = Vec3f(0, -1, 0);
        else if (abs(p.y - c2.y) < 0.01)
            normal = Vec3f(0, 1, 0);
        else if (abs(p.z - c1.z) < 0.01)
            normal = Vec3f(0, 0, -1);
        else if (abs(p.z - c2.z) < 0.01)
            normal = Vec3f(0, 0, 1);

        return normal;
    }
};

struct Plane : Object
{
    Vec3f a, b, c, d;

    Plane(const Vec3f &a, const Vec3f &b, const Vec3f &c, const Vec3f &d, const Material &m) : a(a), b(b), c(c), d(d)
    {
        Object::material = m;
    }

    bool inside(const Vec3f &q) const
    {
        Vec3f n = normal(q);

        Vec3f ua = b - a, ub = c - b, uc = d - c, ud = a - d;
        Vec3f va = q - a, vb = q - b, vc = q - c, vd = q - d;

        if ((cross(ua, va) * n > 0) && (cross(ub, vb) * n > 0) && (cross(uc, vc) * n > 0) && (cross(ud, vd) * n > 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t)
    {
        Vec3f n = normal(p);

        Vec3f vdif = a - p;

        float vdotn = d * n;

        if (fabs(vdotn) < 1.e-4)
            return false;

        float t1 = (vdif * n) / vdotn;

        if (fabs(t) < 0.0001)
            return false;

        Vec3f q = p + d * t;

        if (inside(q))
        {
            t = t1;
            return true;
        }
        else
            return false;
    }

    Vec3f normal(const Vec3f &p) const
    {
        Vec3f result = cross(b - a, c - a);
        result.normalize();
        return result;
    }
};

struct Cylinder : Object
{
    Vec3f cen;
    float r;
    float h;

    Cylinder(const Vec3f &cen, const float &r, const float &h, const Material &m) : cen(cen), r(r), h(h)
    {
        Object::material = m;
    }

    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t)
    {
        Vec3f v = cen - p;

        if (v * d < 0)
        {
            return false;
        }
        else
        {
            float a = (d.x * d.x) + (d.z * d.z);
            float b = 2 * ((d.x * (p.x - cen.x)) + (d.z * (p.z - cen.z)));
            float c = ((p.x - cen.x) * (p.x - cen.x)) + ((p.z - cen.z) * (p.z - cen.z)) - (r * r);

            float diskr = b * b - 4 * (a * c);

            /*if (fabs(diskr) < 0.001)
            {
                //cout << "1" << endl;
                return false;
            }*/

            if (diskr < 0.0)
            {
                return false;
            }
            else
            {
                float t1 = (-b + sqrt(diskr)) / (2 * a);
                float t2 = (-b - sqrt(diskr)) / (2 * a);

                float t3 = min(t1, t2);
                float t4 = max(t1, t2);

                if ((p.y + t3 * d.y >= cen.y) && (p.y + t3 * d.y <= cen.y + h))
                {
                    t = t3;
                    return true;
                }
                else if ((p.y + t4 * d.y >= cen.y) && (p.y + t4 * d.y <= cen.y + h))
                {
                    t = t4;
                    return true;
                }
                else
                    return false;
            }
        }
    }

    Vec3f normal(const Vec3f &p) const
    {
        Vec3f n = (p - cen).normalize();
        n.y = 0;
        return n;
    }
};

struct Triangle : Object
{
    Vec3f v0, v1, v2;
    Vec3f n0, n1, n2;
    Vec3f n;

    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Material &m) : v0(v0), v1(v1), v2(v2)
    {
        Object::material = m;
    }

    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Vec3f &n0, const Vec3f &n1, const Vec3f &n2, const Material &m) : v0(v0), v1(v1), v2(v2), n0(n0), n1(n1), n2(n2)
    {
        Object::material = m;
    }

    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t)
    {
        const float EPSILON = 0.0000001;
        bool triangle_hit;

        Vec3f edge1, edge2, h, s, q;
        float a, f, u, v;
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        h = cross(d, edge2);
        a = edge1 * h;
        if (fabs(a) < EPSILON)
            return false;
        f = 1.0 / a;
        s = p - v0;
        u = f * (s * h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = cross(s, edge1);
        v = f * (d * q);
        if (v < 0.0 || u + v > 1.0)
            return false;

        float intersection = f * (edge2 * q);

        if (intersection > EPSILON)
        {
            t = triangle_hit ? min(t, intersection) : intersection;
            triangle_hit = true;
        }

        //https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/shading-normals
        if (SMOOTH_SHADING)
        {
            n = n0 * (1 - u - v) + n1 * u + n2 * v;
            n.normalize();
        }
        else
        {
            n = cross(v1 - v0, v2 - v0);
            n.normalize();
        }

        return triangle_hit;
    }

    Vec3f normal(const Vec3f &p) const
    {
        return n;
    }
};

#define tupleD tuple<double, double, double>

struct TriangleMesh
{
    Vec3f position;
    Vec3f scale;
    Vec3f rotation;
    Material material;

    vector<OBJParser::Face> faces;
    vector<Triangle> triangles;

    Matrix44f rotMatX;
    Matrix44f rotMatY;
    Matrix44f rotMatZ;
    Matrix44f rotationMatrix;

    TriangleMesh(const char *file, const Material &material) : position(Vec3f(0, 0, 0)), scale(Vec3f(1, 1, 1)), rotation(Vec3f(0, 0, 0)), material(material)
    {
        setup_rotation_matrices();

        OBJParser parser(file);

        faces = parser.get_f();

        create_mesh();
    }

    TriangleMesh(const char *file, const Vec3f &position, const Vec3f &scale, const Vec3f &rotation, const Material &material) : position(position), scale(scale), rotation(rotation), material(material)
    {
        setup_rotation_matrices();

        OBJParser parser(file);

        faces = parser.get_f();

        create_mesh();
    }

    void setup_rotation_matrices()
    {
        rotMatX = Matrix44f(1, 0, 0, 0,
                            0, std::cos(rotation.x), -std::sin(rotation.x), 0,
                            0, std::sin(rotation.x), std::cos(rotation.x), 0,
                            0, 0, 0, 1);

        rotMatY = Matrix44f(std::cos(rotation.y), 0, std::sin(rotation.y), 0,
                            0, 1, 0, 0,
                            -std::sin(rotation.y), 0, std::cos(rotation.y), 0,
                            0, 0, 0, 1);

        rotMatZ = Matrix44f(std::cos(rotation.z), -std::sin(rotation.z), 0, 0,
                            std::sin(rotation.z), std::cos(rotation.z), 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1);

        rotationMatrix = rotMatZ * rotMatY * rotMatX;
    }

    void create_mesh()
    {
        cout << faces.size() << endl;

        for (auto &&face : faces)
        {
            vector<tupleD> verts = face.V;

            Vec3f v0(get<0>(verts[0]) * scale.x + position.x, get<1>(verts[0]) * scale.y + position.y, get<2>(verts[0]) * scale.z + position.z);
            Vec3f v1(get<0>(verts[1]) * scale.x + position.x, get<1>(verts[1]) * scale.y + position.y, get<2>(verts[1]) * scale.z + position.z);
            Vec3f v2(get<0>(verts[2]) * scale.x + position.x, get<1>(verts[2]) * scale.y + position.y, get<2>(verts[2]) * scale.z + position.z);

            rotationMatrix.multVecMatrix(v0, v0);
            rotationMatrix.multVecMatrix(v1, v1);
            rotationMatrix.multVecMatrix(v2, v2);

            auto norms = face.VN;

            Vec3f n0 = Vec3f(get<0>(norms[0]), get<1>(norms[0]), get<2>(norms[0]));
            Vec3f n1 = Vec3f(get<0>(norms[1]), get<1>(norms[1]), get<2>(norms[1]));
            Vec3f n2 = Vec3f(get<0>(norms[2]), get<1>(norms[2]), get<2>(norms[2]));

            rotationMatrix.multVecMatrix(n0, n0);
            rotationMatrix.multVecMatrix(n1, n1);
            rotationMatrix.multVecMatrix(n2, n2);

            Triangle triangle(v0, v1, v2, n0, n1, n2, material);

            triangles.push_back(triangle);
        }
    }

    void insert_to_objects(Objects &objs)
    {
        for (Triangle &triangle : triangles)
        {
            objs.push_back(&triangle);
        }
    }
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const Objects &objs, Vec3f &hit, Material &material, Vec3f &N)
{
    // inicijalno, pretp. da je predmet daleko
    float dist = std::numeric_limits<float>::max();
    float obj_dist = dist;

    for (auto obj : objs)
    {

        if (obj->ray_intersect(orig, dir, obj_dist) && obj_dist < dist)
        {
            dist = obj_dist;
            hit = orig + dir * obj_dist;
            N = obj->normal(hit);
            material = obj->material;
        }
    }

    // provjeri je li sjecište predaleko
    return dist < 1000;
}

Vec3f reflect(const Vec3f &dir, const Vec3f &hit_normal, const Material &hit_material)
{
    return dir + (-hit_normal) * hit_material.refraction_index;
}

// funkcija koja vraca udaljenost sjecista pravca i sfere
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Objects &objs, const Lights &lights, const float &depth, EnvironmentMap &envmap)
{
    int maxDepth = 3;
    Vec3f hit_point;
    Vec3f hit_normal; // normala na povrsinu
    Material hit_material;

    if (!scene_intersect(orig, dir, objs, hit_point, hit_material, hit_normal) || depth > maxDepth)
    {
        return envmap.calc(orig, dir);
    }
    else
    {
        float diffuse_light_intensity = 0;
        float specular_light_intensity = 0;
        float mirroring_intensity = 0.3;

        for (auto light : lights)
        {
            Vec3f light_dir = (light.position - hit_point).normalize();
            float light_dist = (light.position - hit_point).norm();

            Vec3f shadow_orig;
            Vec3f shadow_hit_point;
            Vec3f shadow_hit_normal;
            Material shadow_material;

            if (light_dir * hit_normal < 0)
            {
                hit_normal = -hit_normal;
            }

            shadow_orig = hit_point + hit_normal * 0.001;

            if (scene_intersect(shadow_orig, light_dir, objs, shadow_hit_point, shadow_material, shadow_hit_normal) && (shadow_orig - shadow_hit_point).norm() < light_dist)
                continue;

            diffuse_light_intensity += light.intensity * std::max(0.f, light_dir * hit_normal);

            Vec3f view_dir = (orig - hit_point).normalize();
            Vec3f half_vec = (view_dir + light_dir).normalize();
            specular_light_intensity += light.intensity * powf(std::max(0.f, half_vec * hit_normal), hit_material.specular_exponent);
        }

        Vec3f refraction_vec = reflect(dir, hit_normal, hit_material);

        return hit_material.diffuse_color * hit_material.albedo[0] * diffuse_light_intensity + Vec3f(1, 1, 1) * hit_material.albedo[1] * specular_light_intensity + cast_ray(hit_point + hit_normal * 0.01, (dir - (hit_normal * (dir * hit_normal)) * 2.0), objs, lights, depth + 1, envmap) * mirroring_intensity + cast_ray(hit_point + hit_normal * 0.01, refraction_vec, objs, lights, (hit_material.alpha == 1 ? 13 : depth + 1), envmap) * (1 - hit_material.alpha);
    }
}

//https://www.rose-hulman.edu/class/cs/csse451/examples/notes/present6.pdf
//https://www.scratchapixel.com/lessons/3d-basic-rendering/rendering-3d-scene-overview/perspective-projection
//https://www.scratchapixel.com/lessons/3d-basic-rendering/rendering-3d-scene-overview/visibility-problem
//https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/building-basic-perspective-projection-matrix
Vec3f fromCamToPixel(Vec3f dir, Camera &camera, Viewport &view, int i, int j)
{
    float scale = tan(view.fov / 2.0);
    float distance = (view.nx * 0.5) / scale;

    Vec3f res_dir = camera.dir * distance + camera.u * (i - view.ny * 0.5) + camera.v * (j - view.nx * 0.5);
    res_dir.normalize();

    //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    res_dir = res_dir * cos(camera.rotation_angle) + cross(camera.dir, res_dir) * sin(camera.rotation_angle) + camera.dir * (camera.dir * res_dir) * (1 - cos(camera.rotation_angle));
    res_dir.normalize();

    return res_dir;
}

// funkcija za crtanje
void render(const Objects &objects, const Lights &lights, Viewport &view, Camera &camera, EnvironmentMap &envmap)
{
    vector<Vec3f> buffer(view.nx * view.ny);

    for (int i = 0; i < view.ny; i++)
    {
        for (int j = 0; j < view.nx; j++)
        {
            Vec3f dir = fromCamToPixel(dir, camera, view, i, j);

            buffer[i * view.nx + j] = (cast_ray(camera.pos, dir, objects, lights, 0, envmap));
        }
    }

    ofstream ofs;
    ofs.open("./scena.ppm", ofstream::binary);
    ofs << "P6\n"
        << view.nx << " " << view.ny << "\n255\n";

    for (int i = 0; i < view.ny * view.nx; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            unsigned char color = (unsigned char)(255.f * max(0.f, min(1.f, buffer[i][j])));
            ofs << color;
        }
    }

    ofs.close();
}

int main()
{
    Material red = Material(Vec2f(0.6, 0.3), Vec3f(1, 0, 0), 60, 0.05, 0.7);
    Material green = Material(Vec2f(0.6, 0.3), Vec3f(0, 0.5, 0), 60, 1, 1);
    Material blue = Material(Vec2f(0.9, 0.1), Vec3f(0, 0, 1), 60, 1, 1);
    Material gray = Material(Vec2f(0.9, 0.1), Vec3f(0.5, 0.5, 0.5), 10, 1, 1);
    Material black = Material(Vec2f(0.6, 0.3), Vec3f(0, 0, 0), 60, 1, 1);
    Material yellow = Material(Vec2f(0.6, 0.3), Vec3f(0.94, 0.92, 0.16), 60, 1, 1);
    Material brown = Material(Vec2f(0.6, 0.3), Vec3f(0.29, 0.16, 0.15), 60, 1, 1);

    Objects objs;

    Cuboid surface(Vec3f(-25, -6, -30), Vec3f(25, -5, -9), gray);

    Sphere sphere1(Vec3f(-5, -3.5, -12), 1.5, green);
    Sphere sphere2(Vec3f(3, -4.5, -11.5), 0.5, blue);
    Cuboid cuboid1(Vec3f(1, -5, -15), Vec3f(3, -3, -13), red);

    TriangleMesh djiraf("./objs/djiraf.obj", Vec3f(-14, -5, -11), Vec3f(0.02, 0.02, 0.02), Vec3f(0, 0.785398163, 0), yellow);
    TriangleMesh cilindar("./objs/cilindar.obj", Vec3f(6, -4, -11), Vec3f(1, 1, 1), Vec3f(0, 0, 0), yellow);
    TriangleMesh cilindar2("./objs/cilindar.obj", Vec3f(15, -4.5, -15), Vec3f(1, 0.5, 4), Vec3f(0, 0, 0), green);
    TriangleMesh kocka("./objs/cube.obj", Vec3f(0, -4, -20), Vec3f(1, 1, 1), Vec3f(0, pi / 6, 0), blue);
    TriangleMesh pas("./objs/pas.obj", Vec3f(-12, -5, -15), Vec3f(1, 1, 1), Vec3f(0, 0, 0), brown);

    objs.push_back(&surface);
    objs.push_back(&sphere1);
    objs.push_back(&sphere2);
    objs.push_back(&cuboid1);

    cilindar.insert_to_objects(objs);
    cilindar2.insert_to_objects(objs);
    djiraf.insert_to_objects(objs);
    kocka.insert_to_objects(objs);
    pas.insert_to_objects(objs);

    // definiraj svjetla
    Light l1 = Light(Vec3f(-20, 20, 20), 1);
    Light l2 = Light(Vec3f(20, 30, 20), 1.5);
    Lights lights = {l1, l2};

    Viewport view(1024, 768, pi / 2);
    Camera cam(Vec3f(0, 0, 5), Vec3f(0, 0, -1), 0);

    Viewport view2(700, 700, pi / 4);
    Camera cam2(Vec3f(0, 50, 20), Vec3f(0, -1, -1), 0);

    Viewport view3(768, 1024, 90);
    Camera cam3(Vec3f(0, 0, 0), Vec3f(0, 0, -1), 45);

    Viewport view4(1024, 768, 90);
    Camera cam4(Vec3f(30, 20, 20), Vec3f(-2, 0, -1), 3 * pi / 2);

    EnvironmentMap envmap("./env/envmap2k.ppm", 200);

    render(objs, lights, view, cam, envmap);

    return 0;
}
