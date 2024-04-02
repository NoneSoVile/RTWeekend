#ifndef CAMERA_H
#define CAMERA_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "color.h"
#include "hittable.h"
#include "material.h"

#include <iostream>
#include <png.h>
#include <vector>
#include <jpeglib.h>

class camera {
  public:
    double aspect_ratio      = 1.0;  // Ratio of image width over height
    int    image_width       = 100;  // Rendered image width in pixel count
    int    samples_per_pixel = 10;   // Count of random samples for each pixel
    int    max_depth         = 10;   // Maximum number of ray bounces into scene
    color  background;               // Scene background color

    double vfov     = 90;              // Vertical view angle (field of view)
    point3 lookfrom = point3(0,0,-1);  // Point camera is looking from
    point3 lookat   = point3(0,0,0);   // Point camera is looking at
    vec3   vup      = vec3(0,1,0);     // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    void render(const hittable& world, const hittable& lights, std::string resultImage="result.png") {
        initialize();

        std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        bool useJPEG = false;
        std::string filename;
        std::string substring = ".jpg";

        size_t found = resultImage.find(substring);
        if (found != std::string::npos) {
            useJPEG = true;
            filename = resultImage;
        }
        else {
            found = resultImage.find(".png");
            if (found != std::string::npos) {
                useJPEG = false;
                filename = resultImage;
            }
            else {
                filename = "result.png";
                useJPEG = false;
            }
            
        }

        int quality = 300;
        unsigned char* buf = NULL;
        if(!useJPEG)
            buf = new unsigned char[image_width * image_height * 4 ];
        else
            buf = new unsigned char[image_width * image_height * 6 ];

        for (int j = 0; j < image_height; ++j) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0,0,0);
                color outColor;
                for (int s_j = 0; s_j < sqrt_spp; ++s_j) {
                    for (int s_i = 0; s_i < sqrt_spp; ++s_i) {
                        ray r = get_ray(i, j, s_i, s_j);
                        pixel_color += ray_color(r, max_depth, world, lights);
                    }
                }
                //write_color(std::cout, pixel_color, samples_per_pixel);
                get_color(outColor, pixel_color, samples_per_pixel);
                if(!useJPEG){
                    buf[j * image_width * 4 + i * 4] = int(outColor.e[0]);
                    buf[j*image_width*4 + i*4 + 1] = int(outColor.e[1]);
                    buf[j*image_width*4 + i*4 + 2] = int(outColor.e[2]);
                    buf[j*image_width*4 + i*4 + 3] = 255;
                }else{
                    buf[j*image_width*3 + i*3] = int(outColor.e[0]);
                    buf[j*image_width*3 + i*3 + 1] = int(outColor.e[1]);
                    buf[j*image_width*3 + i*3 + 2] = int(outColor.e[2]);                  
                }
            }
        }

        std::cout << "\nbegin to save result..." << std::endl;


        if(!useJPEG)
            savePNGFile(filename.c_str(), image_width, image_height, buf);
        else
            write_jpeg_file(filename.c_str(), image_width, image_height, buf, quality);
        delete[] buf;
        std::cout << "saved result to " << filename <<  std::endl;

        std::clog << "\rDone.                 \n";
    }

    void savePNGFile(std::string filename, int w, int h, unsigned char* dataBuffer) {
        int width = w;
        int height = h;
        // Create file pointer
        FILE* file = fopen(filename.c_str(), "wb");
        if (!file) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        // Create PNG structures
        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) {
            std::cerr << "Error creating PNG write struct" << std::endl;
            fclose(file);
            return;
        }

        png_infop info = png_create_info_struct(png);
        if (!info) {
            std::cerr << "Error creating PNG info struct" << std::endl;
            png_destroy_write_struct(&png, nullptr);
            fclose(file);
            return
                ;
        }

        // Set up error handling
        if (setjmp(png_jmpbuf(png))) {
            std::cerr << "Error during PNG creation" << std::endl;
            png_destroy_write_struct(&png, &info);
            fclose(file);
            return;
        }

        // Set PNG file I/O
        png_init_io(png, file);

        // Set image properties
        png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

        // Write PNG header
        png_write_info(png, info);

        std::cout << "writing image data" << "\n";
        // Write image data
        png_bytep* rowPointers = new png_bytep[height];
        for (int y = 0; y < height; ++y) {
            rowPointers[y] = png_bytep(dataBuffer + y * width*4);
        }
        png_write_image(png, rowPointers);
        std::cout << "finished writing image data" << "\n";

        // Write end of PNG
        png_write_end(png, nullptr);

        // Clean up
        delete[] rowPointers;
        png_destroy_write_struct(&png, &info);
        fclose(file);

        std::cout << "PNG file saved successfully: " << filename << std::endl;
    }


    void write_jpeg_file(const char* filename, int width, int height, unsigned char* image_buffer, int quality) {
        struct jpeg_compress_struct cinfo;
        struct jpeg_error_mgr jerr;

        FILE* outfile = fopen(filename, "wb");
        if (!outfile) {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            return;
        }

        cinfo.err = jpeg_std_error(&jerr);
        jpeg_create_compress(&cinfo);
        jpeg_stdio_dest(&cinfo, outfile);

        cinfo.image_width = width;
        cinfo.image_height = height;
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_RGB;
        cinfo.next_scanline = 0;

        jpeg_set_defaults(&cinfo);
        jpeg_set_quality(&cinfo, quality, TRUE);

        jpeg_start_compress(&cinfo, TRUE);

         std::cout << "writing jpeg image data" << "\n";
        JSAMPROW row_pointer;
        while (cinfo.next_scanline < cinfo.image_height) {
            row_pointer = (JSAMPROW)&image_buffer[cinfo.next_scanline *
                width * 3];
            jpeg_write_scanlines(&cinfo, &row_pointer, 1);
        }
         std::cout << "finished writing jpeg image data" << "\n";
        jpeg_finish_compress(&cinfo);
        fclose(outfile);
        jpeg_destroy_compress(&cinfo);
    }

  private:
    int    image_height;    // Rendered image height
    int    sqrt_spp;        // Square root of number of samples per pixel
    double recip_sqrt_spp;  // 1 / sqrt_spp
    point3 center;          // Camera center
    point3 pixel00_loc;     // Location of pixel 0, 0
    vec3   pixel_delta_u;   // Offset to pixel to the right
    vec3   pixel_delta_v;   // Offset to pixel below
    vec3   u, v, w;         // Camera frame basis vectors
    vec3   defocus_disk_u;  // Defocus disk horizontal radius
    vec3   defocus_disk_v;  // Defocus disk vertical radius

    void initialize() {
        image_height = static_cast<int>(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta/2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (static_cast<double>(image_width)/image_height);

        sqrt_spp = static_cast<int>(sqrt(samples_per_pixel));
        recip_sqrt_spp = 1.0 / sqrt_spp;

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        vec3 viewport_u = viewport_width * u;    // Vector across viewport horizontal edge
        vec3 viewport_v = viewport_height * -v;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors to the next pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u/2 - viewport_v/2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j, int s_i, int s_j) const {
        // Get a randomly-sampled camera ray for the pixel at location i,j, originating from
        // the camera defocus disk, and randomly sampled around the pixel location.

        auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
        auto pixel_sample = pixel_center + pixel_sample_square(s_i, s_j);

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;
        auto ray_time = random_double();

        return ray(ray_origin, ray_direction, ray_time);
    }

    vec3 pixel_sample_square(int s_i, int s_j) const {
        // Returns a random point in the square surrounding a pixel at the origin, given
        // the two subpixel indices.
        auto px = -0.5 + recip_sqrt_spp * (s_i + random_double());
        auto py = -0.5 + recip_sqrt_spp * (s_j + random_double());
        return (px * pixel_delta_u) + (py * pixel_delta_v);
    }

    vec3 pixel_sample_disk(double radius) const {
        // Generate a sample from the disk of given radius around a pixel at the origin.
        auto p = radius * random_in_unit_disk();
        return (p[0] * pixel_delta_u) + (p[1] * pixel_delta_v);
    }

    point3 defocus_disk_sample() const {
        // Returns a random point in the camera defocus disk.
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
    }

    color ray_color(const ray& r, int depth, const hittable& world, const hittable& lights)
    const {
        hit_record rec;

        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return color(0,0,0);

        // If the ray hits nothing, return the background color.
        if (!world.hit(r, interval(0.001, infinity), rec))
            return background;

        scatter_record srec;
        color color_from_emission = rec.mat->emitted(r, rec, rec.u, rec.v, rec.p);

        if (!rec.mat->scatter(r, rec, srec))
            return color_from_emission;

        if (srec.skip_pdf) {
            return srec.attenuation * ray_color(srec.skip_pdf_ray, depth-1, world, lights);
        }

        auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
        mixture_pdf p(light_ptr, srec.pdf_ptr);

        ray scattered = ray(rec.p, p.generate(), r.time());
        auto pdf_val = p.value(scattered.direction());

        double scattering_pdf = rec.mat->scattering_pdf(r, rec, scattered);

        color sample_color = ray_color(scattered, depth-1, world, lights);
        color color_from_scatter = (srec.attenuation * scattering_pdf * sample_color) / pdf_val;

        return color_from_emission + color_from_scatter;
    }
};


#endif
