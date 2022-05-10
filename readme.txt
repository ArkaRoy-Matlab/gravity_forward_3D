% ***************************************************************
% *** Help file for running all codes
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************
This is a help file for a description of all Data, Source Code, and Subroutine used for the implementation of our present paper 
'Generalized Gauss-FFT 3D forward gravity modeling for irregular topographic mass having any 3D variable density contrast'  

(Copy files and folders in a single folder and run the codes)

	1. Data Files
		(i) Input:
			a.synthetic_topo_fixed_density_shallower_layer
			b.synthetic_topo_fixed_density_deeper_layer
			c.synthetic_topo_exp_density_shallower_layer
			d.synthetic_topo_exp_density_deeper_layer
			e.synthetic_topo_polynomial_density_shallower_layer
			f.synthetic_topo_polynomial_density_deeper_layer
			g.synthetic_topo_complex_density_shallower_layer
			h.synthetic_topo_complex_density_deeper_layer
			i.real_topo_santos
			j.synthetic_x_fixed_density
			k.synthetic_y_fixed_density
			l.synthetic_x_exp_density
			m.synthetic_y_exp_density
			o.synthetic_x_polynomial_density
			p.synthetic_y_polynomial_density
			q.synthetic_x_complex_density
			r.synthetic_y_complex_density
			s.real_x_santos
			t.real_y_santos
		
	Files (a) and (b) are the gridded synthetic topographic data files due to fixed density distributions for shallower and deeper surfaces respectively.  Files (c) and (d) are the gridded synthetic topographic data files due to exponential density distributions for shallower and deeper surfaces respectively. Files (e) and (f) are the gridded synthetic topographic data files due to polynomial density distributions for shallower and deeper surfaces respectively. Files (g) and (h) are the gridded synthetic topographic data files due to complex density distributions for shallower and deeper surfaces respectively. File (i) is the digitized gridded topography data for real sedimentary basin from Chintalpudi sub basin, India. The real data digitized from Chakravarthi, V., and N. Sundararajan, 2007, 3d gravity inversion of basement relief—a depth-dependent density approach: Geophysics, 72, I23–I32. Files (j) and (k) are one dimensional x and y grids for fixed density topographic surfaces. Files (l) and (m) are one dimensional x and y grids for exponential density topographic surfaces. Files (o) and (p) are one dimensional x and y grids for polynomial density topographic surfaces. Files (q) and (r) are one dimensional x and y grids for complex density topographic surfaces. Files (s) and (t) are one dimensional x and y grids for real topographic surface. 

		(ii) Output:
			a1.gravity_fixed_density_prism
			a2.gravity_fixed_density_analytic
			a3.gravity_fixed_density_quadrature1
			a4.gravity_fixed_density_quadrature2
			b1.gravity_exp_density_layer
			b2.gravity_exp_density_analytic
			b3.gravity_exp_density_quadrature1
			b4.gravity_exp_density_quadrature2
			c1.gravity_polynomial_density_layer
			c2.gravity_polynomial_density_analytic
			c3.gravity_polynomial_density_quadrature1
			c4.gravity_polynomial_density_quadrature2
			d1.gravity_complex_density_layer
			d2.gravity_complex_density_quadrature1
			d3.gravity_complex_density_quadrature1
			e1.gravity_real_santos
			e2.gravity_real_density_quadrature
			f1.x_meshgrid_fixed_density
			f2.y_meshgrid_fixed_density
			g1.x_meshgrid_exp_density
			g2.y_meshgrid_exp_density
			h1.x_meshgrid_polynomial_density
			h2.y_meshgrid_polynomial_density		
			i1.x_meshgrid_complex_density
			i2.y_meshgrid_complex_density	
			j1.x_meshgrid_real_density
			j2.y_meshgrid_real_density	

	Files (a1), (a2), (a3) and (a4) are forward gravity anomalies due to fixed density distribution using prismatic approach, analytic Gauss-FFT approach, quadrature based Gauss-FFT approach and quadrature based standard FFT approach respectively. Files (b1), (b2), (b3) and (b4) are forward gravity anomalies due to exponential density distribution using vertically layered approach, analytic Gauss-FFT approach, quadrature based Gauss-FFT approach and quadrature based standard FFT approach respectively. Files (c1), (c2), (c3) and (c4) are forward gravity anomalies due to polynomial density distribution using vertically layered approach, analytic Gauss-FFT approach, quadrature based Gauss-FFT approach and quadrature based standard FFT approach respectively. Files (d1), (d2) and (d3) are forward gravity anomalies due to complex density distribution using vertically layered approach, quadrature based Gauss-FFT approach and quadrature based standard FFT approach respectively. File (e1) and (e2) are the observed gravity anomalies due to Chintanpudi sub basin digitized from Chakravarthi, V., and N. Sundararajan, 2007, 3d gravity inversion of basement relief—a depth-dependent density approach: Geophysics, 72, I23–I32 and quadrature based Gauss-FFT forward gravity anomalies from inverted basement depth. Files  (f1) and (f2) are x and y horizontal observation grid for fixed density model. Files  (g1) and (g2) are x and y horizontal observation grid for exponential density model. Files  (h1) and (h2) are x and y horizontal observation grid for polynomial density model. Files  (i1) and (i2) are x and y horizontal observation grid for complex density model. Files  (j1) and (j2) are x and y horizontal observation grid for real Chintalpudi sub basin.

	2. Matlab functions:
		(i)   center_grid
		(ii)  find_delta
		(iii) find_nodes
		(iv)  gprism
		(v)   grav_layer_gaussfft
		(vi)  grav_quadrature_fft
		(vii) grav_quadrature_gaussfft
		(viii)sfft_X
		(ix)  sfft_Y
		(x)   sfft2 
		(xi)  sifft_X
		(xii) sifft_Y
		(xiii)sifft2 
		(xiv) lgwt
		(xv)  makecolormap

	(i)  center_grid: Matlab function for creating center grid of any 2D data using mean of adjacent grids. 
	(ii) find_delta: Matlab function for finding best delta for any topographic surface for approximating exponential term in Taylor series used in Gauss-FFT forward Modelling. 
	
	(iii) find_nodes: Matlab function for finding best number of Gaussian node for quadrature integral for given x grid, y grid, gaussian quadrature weight and density distributions. 

	(iv) gprism: Analytical solution to compute the vertical attraction of a rectangular prism. Sides of prism are parallel to x,y,z axes, and z is vertical down. This subroutine is adopted from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68. 
	(v)   grav_layer_gaussfft: Matlab function for finding gravity anomalies for topographic surfaces having any 3D density contrast using vertially sliced horizontal layers.

	(vi) grav_quadrature_fft: Matlab function for finding gravity anomalies for topographic surfaces having any 3D density contrast quadrature based standard FFT method.

	vii) grav_quadrature_gaussfft: Matlab function for finding gravity anomalies for topographic surfaces having any 3D density contrast using quadrature based Gauss-FFT method.

	(viii) sfft_X: 1D shift FFT along each row of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(ix) sfft_Y: 1D shift FFT along each column of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(x) sfft2 : 2D shift FFT of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(xi) sifft_X: 1D shift inverse FFT along each row of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(xii) sifft_Y: 1D shift inverse FFT along each column of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(xiii) sifft2 : 2D shift inverse FFT of any 2D matrix. Inspired from Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68.

	(xiv)  lgwt.m - This script is for computing definite integrals using Legendre-Gauss 
 Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval [a,b] with truncation order N. Suppose you have a continuous function f(x) which is defined on [a,b]
which you can evaluate at any x in [a,b]. Simply evaluate it at all of the values contained in the x vector to obtain a vector f. Then compute the definite integral using sum(f.*w);

	This code is written by Greg von Winckel - 02/25/2004. Here we have used it for our calculation. 
	
	(xv) makecolormap.m - Matlab function for creating colormap as per user choice, all color name and corresponding color representation can be found from 
the web address  http://kodama.fubuki.info/wiki/wiki.cgi/GrADS/script/color.gs?lang=en

	3. Matlab scripts:
		(i)    all_model_layer
		(ii)   all_model_quadrature
		(iii)  prism_fixed_density
		(iv)   fixed_density_analytic
		(v)    exponential_density_analytic
		(vi)   polynomial_density_analytic
		(vii)  plllt1
		(viii) plllt2
		(ix)   generalized

	(i) all_model_layer: Matlab code for finding gravity anomalies for different topographic surfaces and corresponding density contrasts using vertically sliced horizontal layers.
	
	(ii) all_model_quadrature: Matlab code for finding gravity anomalies for different topographic surfaces and corresponding density contrasts using Gauss-FFT and standard FFT quadrature based model

	(iii) prism_fixed_density: Matlab code for prism model for finding gravity anomalies of topographic mass having fixed density distribution using prismatic approach.

	(iv) fixed_density_analytic: Matlab code for finding gravity anomaly due uniform density distribution using analytically derived integral and Gauss-FFT method.

	(v) exponential_density_analytic: Matlab code for finding gravity anomaly due exponential density distribution using analytically derived integral and Gauss-FFT method.

	(vi) polynomial_density_analytic: Matlab code for finding gravity anomaly due polynomial density distribution using analytically derived integral and Gauss-FFT method.

	(vii) plllt1: Matlab code for plotting all synthetic topographic surfaces and corresponding analytic forward gravity anomalies.

	(viii) plllt2: Matlab code for plotting all synthetic topographic surfaces and corresponding gauss FFT and standard FFT quadrature based forward gravity anomalies.

  	(ix) generalized: A generalized matlab code to evaluate gravity anomalies by any user by providing the input data. 	

    by running this code plllt1.m, you can reproduce Figures 2, 4, 6, 8, 10, and 11 in the manuscript.
	by running this code plllt2.m, you can reproduce Figures 3, 5, 7, and 9 in the manuscript.

