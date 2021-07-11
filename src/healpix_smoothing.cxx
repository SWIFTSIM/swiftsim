#include <stdio.h>
#include <math.h>

#include "healpix_cxx/vec3.h"
#include "healpix_cxx/pointing.h"
#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/datatypes.h"

extern "C" {
#include "atomic.h"
#include "projected_kernel.h"
}

// 2D Wendland C2 kernel (omitting normalisation)
// TODO: use same kernel as Swift?
static double projected_kernel(double r, double h) {
  
  const double w = pow(1.0-r/h, 4.0)*(1.0+4*r/h);
  if(w > 0)
    return w;
  else
    return 0.0;
}


struct healpix_smoothing_info {
  int nside;
  double max_pixrad;
  double kernel_gamma;
  Healpix_Base2 healpix_base;
  struct projected_kernel_table kernel;
};


extern "C" {

  struct healpix_smoothing_info *healpix_smoothing_init(int nside, double gamma) {
    
    const Healpix_Ordering_Scheme scheme = RING;
    healpix_smoothing_info *smooth_info = new struct healpix_smoothing_info;
    smooth_info->nside = nside;
    smooth_info->healpix_base = Healpix_Base2(nside, scheme, SET_NSIDE);
    smooth_info->max_pixrad = smooth_info->healpix_base.max_pixrad();
    smooth_info->kernel_gamma = gamma;
    projected_kernel_init(&smooth_info->kernel);
    return smooth_info;
  }

  void healpix_smoothing_clean(struct healpix_smoothing_info *smooth_info) {
    projected_kernel_clean(&smooth_info->kernel);
    delete smooth_info;
  }

  size_t healpix_smoothing_get_npix(struct healpix_smoothing_info *smooth_info) {
    return (size_t) smooth_info->healpix_base.Npix();
  }

  double healpix_smoothing_get_max_pixrad(struct healpix_smoothing_info *smooth_info) {
    return smooth_info->healpix_base.max_pixrad();
  }

  size_t healpix_smoothing_vec2pix(struct healpix_smoothing_info *smooth_info, 
                                   const double *pos) {
    
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    return (size_t) smooth_info->healpix_base.vec2pix(part_vec);
  }

  void healpix_smoothing_vec2ang(struct healpix_smoothing_info *smooth_info,
                                 const double *pos, double *theta, double *phi) {
    
    pointing p = pointing(vec3(pos[0], pos[1], pos[2]));
    *theta = p.theta;
    *phi = p.phi;
  }

  size_t healpix_smoothing_ang2pix(struct healpix_smoothing_info *smooth_info,
                                   const double theta, const double phi) {
    pointing p = pointing(theta, phi);
    return (size_t) smooth_info->healpix_base.ang2pix(p);
  }

  void healpix_smoothing_get_pixel_range(struct healpix_smoothing_info *smooth_info,
                                         const double theta, const double phi, 
                                         const double radius, size_t *first_pixel,
                                         size_t *last_pixel) {
    pointing p = pointing(theta, phi);

    // Input radius is the angle corresponding to the smoothing length,
    // so we need to find all pixels within kernel_gamma*radius
    const double search_radius = radius*smooth_info->kernel_gamma;

    // Small particles get added to a single pixel
    if(search_radius < smooth_info->max_pixrad) {
      int64 pixel = smooth_info->healpix_base.ang2pix(p);
      *first_pixel = pixel;
      *last_pixel = pixel;
      return;
    }

    // Find all pixels with centres within the angular radius
    std::vector<int64> pixels;
    smooth_info->healpix_base.query_disc(p, search_radius, pixels);
    *first_pixel = pixels[0];
    *last_pixel = pixels[0];
    for(size_t i=0; i < pixels.size(); i++) {
      if(pixels[i] < *first_pixel)*first_pixel = pixels[i];
      if(pixels[i] > *last_pixel)*last_pixel = pixels[i];
    }

    return;
  }

  void healpix_smoothing_add_to_map(struct healpix_smoothing_info *smooth_info,
                                    const double theta, const double phi, const double radius,
                                    const double value, const size_t local_pix_offset,
                                    const size_t local_nr_pix, double *map_data) {

    pointing p = pointing(theta, phi);

    // Input radius is the angle corresponding to the smoothing length,
    // so we need to find all pixels within kernel_gamma*radius
    const double search_radius = radius*smooth_info->kernel_gamma;
  
    // Small particles get added to a single pixel
    if(search_radius < smooth_info->max_pixrad) {
      int64 pixel = smooth_info->healpix_base.ang2pix(p);
      if((pixel >= local_pix_offset) && (pixel < local_pix_offset+local_nr_pix))
        atomic_add_d(&map_data[pixel-local_pix_offset], value);
      return;
    }

    // Find all pixels with centres within the angular radius
    std::vector<int64> pixels;
    smooth_info->healpix_base.query_disc(p, search_radius, pixels);
  
    // Check for the case where a particle was sent to an MPI rank it doesn't contribute to
    const size_t npix = pixels.size();
    if(npix == 0)return;

    // Vector to store pixel weights
    std::vector<double> weight(npix);

    // Particle direction vector
    vec3 part_vec = p.to_vec3();

    // Loop over pixels within the radius
    double tot = 0.0;
    for(size_t i=0; i<npix; i+=1) {

      // Get direction vector to centre of this pixel
      vec3 pixel_vec = smooth_info->healpix_base.pix2vec(pixels[i]);

      // Find angle between this pixel centre and the particle
      const double angle = acos(dotprod(pixel_vec, part_vec));

      // Evaluate the kernel at this radius
      weight[i] = projected_kernel_eval(&smooth_info->kernel, angle/radius);
      tot += weight[i];

    }
  
    // Now accumulate contributions to pixels
    for(size_t i=0; i<npix; i+=1) {
      if((pixels[i] >= local_pix_offset) && (pixels[i] < local_pix_offset+local_nr_pix)) {
        atomic_add_d(&(map_data[pixels[i]-local_pix_offset]), weight[i]/tot*value);
      }
    }
  }

} // extern "C"



