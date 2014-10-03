/* W B Langdon at MUN 10 May 2007
 * Program to demonstarte use of OpenGL's glDrawPixels
 */

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include <iostream>
#include <sstream>
#include "math.h"

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

#include <thread>
#include <mutex>


using std::cerr;
using std::cout;
using std::endl;

using boost::variate_generator;
using boost::mt19937;
using boost::exponential_distribution;
using boost::gamma_distribution;
using boost::uniform_real;

unsigned int window_width, window_height;
int window_size;
int num_active_pixels;
int num_threads;

//#define FAST
//#define THREADING
//#define NUMBER_DROPLETS 50

bool pixel_frozen = false;

float *pixels;

size_t *active_pixels;
bool *frozen_pixels;

#ifdef THREADING
std::vector<std::thread*> threads;
std::mutex rng_mtx;
std::mutex freeze_mtx;
#endif

variate_generator<mt19937, uniform_real<> > random_0_1( mt19937(time(0)), uniform_real<>(0.0, 1.0));

/**
 *  Determine the position of a pixel
 */
static int POSITION(int x, int y) {
    return (y * window_width) + x;
}

/**
 *  This will set the given pixel to black
 */
void unset_pixel(int x, int y) {
    pixels[ POSITION(x, y) * 3]         = 0.0;
    pixels[(POSITION(x, y) * 3) + 1]    = 0.0;
    pixels[(POSITION(x, y) * 3) + 2]    = 0.0;
}


/**
 *  This will set the given pixel to white
 */
void set_pixel(int x, int y) {
    pixels[ POSITION(x, y) * 3]         = 1.0;
    pixels[(POSITION(x, y) * 3) + 1]    = 1.0;
    pixels[(POSITION(x, y) * 3) + 2]    = 1.0;
}

/**
 *  This will freeze the given pixel
 */
void set_frozen(int x, int y) {
    frozen_pixels[POSITION(x, y)] = true;
}

/**
 * This checks if the given coordinates are adjacent to a frozen pixel
 *
 *
 */
bool is_adjacent(int x, int y) {
    if (
            (x < window_width - 1 && frozen_pixels[POSITION(x+1, y)]) ||
            (x > 0 && frozen_pixels[POSITION(x-1, y)]) ||
            (y < window_height - 1 && frozen_pixels[POSITION(x, y+1)]) ||
            (y > 0 && frozen_pixels[POSITION(x, y-1)])
            ) {
        return true;
    }
    return false;
}

/**
 * This sets the given x and y values to a random point on the perimeter
 *
 */
void set_start_pos(size_t *x, size_t *y) {
    size_t x_val = 0;
    size_t y_val = 0;

#ifdef THREADING
    rng_mtx.lock();
#endif
    double rand = random_0_1();
#ifdef THREADING
    rng_mtx.unlock();
#endif

    if (rand < static_cast<double>(window_width)/(window_width + window_height)) {
#ifdef THREADING
        rng_mtx.lock();
#endif
        x_val = random_0_1() * window_width;
        rand = random_0_1();
#ifdef THREADING
        rng_mtx.unlock();
#endif
        if (rand < 0.5) {
            y_val = 0;
        } else {
            y_val = window_height - 1;
        }
    } else {
#ifdef THREADING
        rng_mtx.lock();
#endif
        y_val = random_0_1() * window_height;
        rand = random_0_1();
#ifdef THREADING
        rng_mtx.unlock();
#endif
        if (rand < 0.5) {
            x_val = 0;
        } else {
            x_val = window_width - 1;
        }
    }
    *x = x_val;
    *y = y_val;
}

/**
 * This moves every pixel in one of four directions and check if it is next to
 * a frozen pixel
 */
void move_pixels(size_t start, size_t end) {
#ifdef THREADING
    while(true) {
#endif
        for (int i = start; i < end; i++) {
#ifdef THREADING
            rng_mtx.lock();
#endif
            double rand_val = random_0_1();
#ifdef THREADING
            rng_mtx.unlock();
#endif
            int x_pos = i*2;
            int y_pos = i*2 + 1;
            unset_pixel(active_pixels[x_pos], active_pixels[y_pos]);
            if (rand_val < 0.25) { // Move up
                if (active_pixels[x_pos] < window_width - 1) {
                    active_pixels[x_pos] += 1;
                } else {
                    active_pixels[x_pos] -= 1;
                }
            } else if (rand_val < 0.50) { // Move down
                if (active_pixels[x_pos] > 0) {
                    active_pixels[x_pos] -= 1;
                } else {
                    active_pixels[x_pos] += 1;
                }
            } else if (rand_val < 0.75) { // Move left
                if (active_pixels[y_pos] < window_height - 1) {
                    active_pixels[y_pos] += 1;
                } else {
                    active_pixels[y_pos] -= 1;
                }
            } else { // Move right
                if (active_pixels[y_pos] > 0) {
                    active_pixels[y_pos] -= 1;
                } else {
                    active_pixels[y_pos] += 1;
                }
            }
            set_pixel(active_pixels[x_pos], active_pixels[y_pos]);
#ifdef THREADING
            freeze_mtx.lock();
#endif
            if (is_adjacent(active_pixels[x_pos], active_pixels[y_pos])) {
                set_frozen(active_pixels[x_pos], active_pixels[y_pos]);
                set_start_pos(&active_pixels[x_pos], &active_pixels[y_pos]);
                pixel_frozen = true;
            }
#ifdef THREADING
            freeze_mtx.unlock();
#endif
        }
#ifdef THREADING
    }
#endif
}

void idle() {
#ifndef THREADING
    move_pixels(0, num_active_pixels);
#endif

#ifdef FAST
    if (pixel_frozen) {
#endif
        glutPostRedisplay();
        pixel_frozen = false;
#ifdef FAST
    }
#endif
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //glDrawPixels writes a block of pixels to the framebuffer.
    glDrawPixels(window_width, window_height, GL_RGB, GL_FLOAT, pixels);

    glFlush();
    glutSwapBuffers();
}

int main(int argc, char** argv) {
    if (argc != 4) {
        cerr << "Error running " << argv[0] << ", wrong number of arguments." << endl;
        cerr << "Usage:" << endl;
        cerr << "\t" << "./" << argv[0] << " <window width> <window height> <num particles>" << endl;

        exit(0);
    }

    /**
     *  Set the global variables for the window width, height and size
     */
    window_width = atoi(argv[1]);
    window_height = atoi(argv[2]);
    window_size = window_width * window_height;

    num_active_pixels = atoi(argv[3]);
#ifdef THREADING
    num_threads = std::thread::hardware_concurrency() - 1; // Save one for the main loop

    if (num_active_pixels < num_threads) {
        num_threads = num_active_pixels;
    }
    num_threads = 1;
#endif

    /**
     *  Initialize the pixel matrix
     */
    pixels = (float*)malloc(sizeof(float) * window_size * 3);       //pixles are R G B, so 3 values per pixel
    memset(pixels, 0, window_size * 3);                             // Set every pixel value to zero;
    //for (int i = 0; i < window_size; i++) pixels[i] = 0.0;        //setting R, G and B to 0 will make every pixel black

    frozen_pixels = (bool*)malloc(sizeof(bool) * window_size);      // Create boolean val for every pixel
    memset(frozen_pixels, 0, window_size);                          // set all boolean values to false

    active_pixels = (size_t*)malloc(sizeof(size_t) * num_active_pixels * 2);
    for (int i = 0; i < num_active_pixels; i++) {
        set_start_pos(&active_pixels[i*2], &active_pixels[i*2 + 1]);
        set_pixel(active_pixels[i*2], active_pixels[i*2 + 1]);
    }

    // Set middle pixel to white
    // TODO Why is this backwards?
    set_pixel(window_width/2, window_height/2);
    set_frozen(window_width/2, window_height/2);

    cout << "Initialized snowflake matrix!" << endl;
    cout << "window width: "    << window_width << endl;
    cout << "window height: "   << window_height << endl;
    cout << "window size : "    << window_size << endl;
    cout << "active pixels: "   << num_active_pixels << endl;
#ifdef THREADING
    cout << "num threads: "     << num_threads << " (plus one for main loop)" << endl;
#endif

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("Snowflake Simulation");

    glutDisplayFunc(display);
    //glutReshapeFunc(reshape);
    //glutMouseFunc(mouse_button);
    //glutMotionFunc(mouse_motion);
    //glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    //glPointSize(2);

    //Start threads to move pixels
#ifdef THREADING
    for (int i = 0; i < num_threads; i++) {
        threads.push_back(new std::thread(move_pixels, num_active_pixels*(i/num_threads), num_active_pixels*((i+1)/num_threads)));
    }
#endif

    glutMainLoop();
}
