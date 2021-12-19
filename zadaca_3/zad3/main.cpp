#include <tuple>
#include "geometry.h"
#include "tgaimage.h"
#include "OBJParser.h"

#define pi 3.14159265359

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);

const uint32_t Width = 640;
const uint32_t Height = 480;

// pogled sa sredine
/*const Matrix44f mCam = {1, 0, 0, 0, 
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, -20, 1};*/

// malo rotirani pogled
const Matrix44f mCam = {1, -1, 0, 0,
                        0, 1, 0, 0,
                        -1, -1, 1, 0,
                        0, -1, -20, 1};

const float nearPlane = 1;
const float farPLane = 1000;

const float FOV = 45;

void convertToRaster(const Vec3f &po, const Matrix44f &mCam, const float &l, const float &r, const float &t, const float &b, Vec3f &ps)
{
    Vec3f pCam;

    mCam.multVecMatrix(po, pCam);

    // pretvori u prostor ekrana
    Vec2f pScreen;
    pScreen.x = nearPlane * pCam.x / -pCam.z;
    pScreen.y = nearPlane * pCam.y / -pCam.z;

    // pretvori u kanonski volumen [-1, 1]
    Vec2f pCan;
    pCan.x = 2 * pScreen.x / (r - l) - (r + l) / (r - l);
    pCan.y = 2 * pScreen.y / (t - b) - (t + b) / (t - b);

    // pretvori u raster prostor
    ps.x = (pCan.x + 1) / 2 * Width;

    // invertaj y (u raster prostoru je y prema dole)
    ps.y = (1 - pCan.y) / 2 * Height;
    ps.z = -pCam.z;
}

float edgeFunction(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{
    return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

void triangle(Vec3f &v0, Vec3f &v1, Vec3f &v2, TGAImage &image, TGAImage &texture, float *depthBuffer)
{
    v0.z = 1 / v0.z,
    v1.z = 1 / v1.z,
    v2.z = 1 / v2.z;

    float xmin = min(v0.x, min(v1.x, v2.x));
    float ymin = min(v0.y, min(v1.y, v2.y));
    float xmax = max(v0.x, max(v1.x, v2.x));
    float ymax = max(v0.y, max(v1.y, v2.y));

    // trokut izvan ekrana
    if (xmin > Width - 1 || xmax < 0 || ymin > Height - 1 || ymax < 0)
        return;

    unsigned int x0 = max(0, (int)floor(xmin));
    unsigned int x1 = min(int(Width) - 1, (int)floor(xmax));
    unsigned int y0 = max(0, (int)floor(ymin));
    unsigned int y1 = min(int(Height) - 1, (int)floor(ymax));

    float area = edgeFunction(v0, v1, v2);

    for (uint32_t y = y0; y <= y1; ++y)
    {
        for (uint32_t x = x0; x <= x1; ++x)
        {
            Vec3f pixel(x + 0.5, y + 0.5, 0);

            float w0 = edgeFunction(v1, v2, pixel);
            float w1 = edgeFunction(v2, v0, pixel);
            float w2 = edgeFunction(v0, v1, pixel);

            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                w0 /= area;
                w1 /= area;
                w2 /= area;

                float t = v0.z * w0 + v1.z * w1 + v2.z * w2;
                float z = 1 / t;

                // test dubine
                if (z < depthBuffer[y * Width + x])
                {
                    depthBuffer[y * Width + x] = z;

                    Vec3f v0Cam, v1Cam, v2Cam;
                    mCam.multVecMatrix(v0, v0Cam);
                    mCam.multVecMatrix(v1, v1Cam);
                    mCam.multVecMatrix(v2, v2Cam);

                    float px = (v0Cam.x / -v0Cam.z) * w0 + (v1Cam.x / -v1Cam.z) * w1 + (v2Cam.x / -v2Cam.z) * w2;
                    float py = (v0Cam.y / -v0Cam.z) * w0 + (v1Cam.y / -v1Cam.z) * w1 + (v2Cam.y / -v2Cam.z) * w2;

                    Vec3f pt(px * z, py * z, -z);

                    Vec3f n = (v1Cam - v0Cam).crossProduct(v2Cam - v0Cam);
                    n.normalize();
                    Vec3f dir = -pt;
                    dir.normalize();

                    int u = round(texture.get_width() * float(x - xmin) / (xmax - xmin));
                    int v = round(texture.get_height() * float(y - ymin) / (ymax - ymin));
                    TGAColor texColor = texture.get(u, v);

                    image.set(x, y, texColor);
                }
            }
        }
    }
}

void triangleMesh(string fileName, Vec3f position, Vec3f scale, float l, float r, float t, float b, TGAImage &image, TGAImage &texture, float *depthBuffer)
{
    OBJParser parser(fileName);

    auto faces = parser.get_f();

    for (auto face : faces)
    {
        auto vertices = face.V;

        Vec3f v0 = Vec3f(get<0>(vertices[0]) * scale.x + position.x, get<1>(vertices[0]) * scale.y + position.y, get<2>(vertices[0]) * scale.z + position.z);
        Vec3f v1 = Vec3f(get<0>(vertices[1]) * scale.x + position.x, get<1>(vertices[1]) * scale.y + position.y, get<2>(vertices[1]) * scale.z + position.z);
        Vec3f v2 = Vec3f(get<0>(vertices[2]) * scale.x + position.x, get<1>(vertices[2]) * scale.y + position.y, get<2>(vertices[2]) * scale.z + position.z);

        Vec3f v0Raster, v1Raster, v2Raster;
        convertToRaster(v0, mCam, l, r, t, b, v0Raster);
        convertToRaster(v1, mCam, l, r, t, b, v1Raster);
        convertToRaster(v2, mCam, l, r, t, b, v2Raster);

        triangle(v0Raster, v1Raster, v2Raster, image, texture, depthBuffer);
    }
}

int main()
{
    // definiraj sliku
    TGAImage image(Width, Height, TGAImage::RGB);

    TGAImage texture;
    texture.read_tga_file("./textures/tekstura.tga");

    TGAImage djirafTex;
    djirafTex.read_tga_file("./textures/djiraf.tga");

    TGAImage wallTex;
    wallTex.read_tga_file("./textures/wall.tga");

    Matrix44f mCam = mCam.inverse();

    float left, right, top, bottom;

    float scale = tan(FOV / 2.);

    top = scale * nearPlane;
    right = scale * nearPlane;

    bottom = -top;
    left = -right;

    float *depthBuffer = new float[Width * Height];

    for (int i = 0; i < Width * Height; ++i)
    {
        depthBuffer[i] = farPLane;
    }

    triangleMesh("./objs/pas.obj", Vec3f(-1, -2, 1), Vec3f(1, 1, 1), left, right, top, bottom, image, texture, depthBuffer);
    triangleMesh("./objs/djiraf.obj", Vec3f(2, -2, 0), Vec3f(0.02, 0.02, 0.02), left, right, top, bottom, image, djirafTex, depthBuffer);
    triangleMesh("./objs/cube.obj", Vec3f(-4, 0, 0), Vec3f(1, 1, 1), left, right, top, bottom, image, wallTex, depthBuffer);

    // spremi sliku
    //image.flip_vertically();
    image.write_tga_file("scena.tga");

    delete[] depthBuffer;

    return 0;
}