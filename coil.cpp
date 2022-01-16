#include <cmath>
#include <cstdio>
#include <vector>
#include <utility>
using namespace std;

#include <Magick++.h>
using namespace Magick;

#include <cairo.h>

using pdd = pair<double, double>;

const double PI = acos(-1);

// size of input image
int IMG_WIDTH;
int IMG_HEIGHT;
// size of result
int DST_WIDTH  = 2000;
int DST_HEIGHT = 2000;

// Starting point of the spiral.  (0.5, 0.5) is the center of the canvas.
// Other coordinates may apply, e.g. (1.5, 1.5), (-0.5, 0.5)
// The visible portion of the canvas is [0, 1] x [0, 1].
const double CENTER_X = 0.5;
const double CENTER_Y = 0.5;

const int loop_num = 60; // number of spirals originating from the center
const int seg_base = 40; // divide the first spiral into seg_base pieces
// The i-th spiral gets segmented into (seg_base + i * seg_inc) pieces.
const double seg_inc = 2;

const pdd CORNERS[] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}}; // four corners of the visible canvas
// Equation for the archimedean spiral: phi = a * theta
double a = 0.5 / (2 * PI * loop_num);
// const double a = max(max(CENTER_X, 1 - CENTER_X),
//                      max(CENTER_Y, 1 - CENTER_Y)) / (2 * PI * loop_num);

// Distance between point a & b.
double dist(pdd a, pdd b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return sqrt(dx*dx + dy*dy);
}

// Convert rgb to grayscale (0 means black).
double rgb_to_grayscale(double r, double g, double b) {
    return 0.299 * r + 0.587 * g + 0.114 * b;
}

// Read image into dst_img.
// dst_img is a vector of doubles (1 means black).
void read_image(const char *filename, vector<double> &dst_img) {
    Image img(filename);
    // Magick++ uses either 8 or 16 bit for each r g b.
    const int BIT_DEPTH = sizeof(PixelPacket::red) * 8;
    // either 255 or 65535
    const int PIXEL_MAX = (1 << BIT_DEPTH) - 1;

    IMG_WIDTH = img.columns();
    IMG_HEIGHT = img.rows();

    printf("Original image %s, width: %d, height: %d.\n",
           filename, IMG_WIDTH, IMG_HEIGHT);
    printf("Magick++ color depth: %u bits.\n", BIT_DEPTH);

    // read raw pixels of the image
    const PixelPacket *pixels = img.getConstPixels(0, 0, IMG_WIDTH, IMG_HEIGHT);
    printf("(0, 0): %u %u %u\n", pixels[0].red, pixels[0].green, pixels[0].blue);

    // convert raw rgb pixels to inverted grayscale (1 means black)
    dst_img.resize(IMG_WIDTH * IMG_HEIGHT);
    auto it = dst_img.begin();
    for (int x = 0; x < IMG_HEIGHT; ++x) {
        for (int y = 0; y < IMG_WIDTH; ++y) {
            double grayscale = rgb_to_grayscale(
                pixels->red, pixels->green, pixels->blue) / PIXEL_MAX;
            *it++ = 1 - grayscale; // 1 means black
            ++pixels;
        }
    }
}

// Calculate d1 & d2 based on how dark the img is.
// The region of interest is indicated by
//   r: radius
//   d: with of a single side
//   t0, dt: degree range, [t0, t0 + dt].
// The image is scaled to its longest side (IMG_DIMENSION).
pdd calc_dimensions(const vector<double> &img,
                    double r,
                    double d, // width of a single side (pi * a)
                    double t0,
                    double dt) {
    const int IMG_DIMENSION = max(IMG_WIDTH, IMG_HEIGHT);
    cairo_surface_t *surface = cairo_image_surface_create(
        CAIRO_FORMAT_RGB24, IMG_DIMENSION, IMG_DIMENSION);
    cairo_t *cr = cairo_create(surface);

    cairo_scale(cr, IMG_DIMENSION, IMG_DIMENSION); // square

    // background is black, arc is white
    // the non-zero pixels are part of the arc
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_arc(cr, CENTER_X, CENTER_Y, r, t0, t0 + dt);
    cairo_set_line_width(cr, d * 2);
    cairo_stroke(cr);

    cairo_surface_flush(surface);

    // How many bytes each row takes.
    int stride = cairo_image_surface_get_stride(surface);
    // How many bytes each pixel takes. (I got 4)
    int increment = stride / cairo_image_surface_get_width(surface);

    double area = 0;            // area of the white arc
    double sum = 0;             // arc * image

    unsigned char * const pimg = cairo_image_surface_get_data(surface);
    auto it = img.cbegin();

    // used to center image
    int x0 = (IMG_DIMENSION - IMG_HEIGHT) / 2;
    int y0 = (IMG_DIMENSION - IMG_WIDTH) / 2;

    for (int x = 0; x < IMG_HEIGHT; ++x) {
        for (int y = 0; y < IMG_WIDTH; ++y) {
            // pointer to pixel of size `increment'
            auto p = pimg + (x0 + x) * stride + (y0 + y) * increment;
            unsigned char r = p[0];
            unsigned char g = p[1];
            unsigned char b = p[2];

            // since arc is white, the resulting `grayscale' is right (1 means black)
            double grayscale = rgb_to_grayscale(r, g, b) / 255;
            area += grayscale;
            sum += (*it++) * grayscale;

            // if (r || g || b) {
            //     printf("(%d, %d): %.2lf\n", x, y, grayscale);
            // }
        }
    }
    // For simplicity, make d1 = d2.
    double ratio = max(0.02, sum / max(area, 1.0));
    double d1 = d * ratio;
    double d2 = d * ratio;

    // printf("sum: %.2lf, area: %.2lf, ratio: %.2lf\n", sum, area, ratio);
    // printf("  d1: %e, d2: %e\n", d1, d2);

    // cairo_surface_write_to_png(surface, "arc.png");
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    return {d1, d2};
}

// Convert polar coordinate to cartesian coordinate (used in cairo).
pdd get_point(double th, double r) {
    return {r * cos(th) + CENTER_X, r * sin(th) + CENTER_Y};
}

pair<pdd, pdd> arc(double r, double d1, double d2, double t0, double dt) {
    double r1 = r - d1;
    double r2 = r + d2;
    double th = t0 + dt / 2;
    return {get_point(th, r1), get_point(th, r2)};
}

void draw(cairo_t *cr, const vector<double> &img) {
    // white background
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    // black paint
    cairo_set_source_rgb(cr, 0, 0, 0);

    vector<pdd> p1, p2;

    printf("Drawing archimedean spirals.\n");
    printf("  loop num: %d\n", loop_num);
    printf("  seg base: %d\n", seg_base);
    printf("  seg inc:  %.2lf\n", seg_inc);
    printf("  a:        %e\n", a);
    printf("\n");

    for (int loop = 0; loop < loop_num; ++loop) {
        putchar('*'); fflush(stdout);

        int seg_num = seg_base + loop * seg_inc;
        double dt = 2 * PI / seg_num;
        for (int i = 0; i < seg_num; ++i) {
            double th = 2 * PI * (loop + 1.0 * i / seg_num);
            double r = a * th;
            auto dpair = calc_dimensions(img, r, PI*a, th, dt);
            const double &d1 = dpair.first;
            const double &d2 = dpair.second;
            // printf("loop: %d, i: %d, d1: %lf, d2: %lf\n", loop, i, d1, d2);
            auto ppair = arc(r, d1, d2, th, dt);
            // each portion provides 2 points
            p1.push_back(ppair.first);
            p2.push_back(ppair.second);
        }
    }
    putchar('\n');

    // connect all the points
    cairo_move_to(cr, CENTER_X, CENTER_Y);
    for (auto it = p1.cbegin(); it != p1.cend(); ++it) {
        cairo_line_to(cr, it->first, it->second);
    }
    for (auto it = p2.crbegin(); it != p2.crend(); ++it) {
        cairo_line_to(cr, it->first, it->second);
    }
    cairo_close_path(cr);
    cairo_fill(cr);
}

int main(int argc, char *argv[]) {
    InitializeMagick(nullptr);

    // This makes the coils fill the canvas.
    // for (auto p : CORNERS) {
    //     a = max(a, dist({CENTER_X, CENTER_Y}, p) / (2 * PI * loop_num));
    // }

    // a = 9.378295e-04;           // for 9pics with loop num 360

    const char *filename = argc > 1 ? argv[1] : "img.jpg";

    // calc_dimensions(img, 0.3, 0.01, -PI/4, PI/8);

    vector<double> img;
    read_image(filename, img);

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24,
                                                          DST_WIDTH,
                                                          DST_HEIGHT);
    cairo_t *cr = cairo_create(surface);
    cairo_scale(cr, DST_WIDTH, DST_HEIGHT);

    draw(cr, img);

    cairo_surface_write_to_png(surface, "result.png");
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}
