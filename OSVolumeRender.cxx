#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkImageData.h"
#include "vtkMetaImageReader.h"
#include "vtkCamera.h"
#include "vtkVolume.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

#define WSIZE 1024

int main( int argc, char * argv[] )
{

	char *  in_file_name = argv[1];
	char * out_file_name = argv[2];
	double vux = atof(     argv[3] );
	double vuy = atof(     argv[4] );
	double vuz = atof(     argv[5] );
	double px  = atof(     argv[6] );
	double py  = atof(     argv[7] );
	double pz  = atof(     argv[8] );
  double min;
  double max;
  if ( argc > 9 )
    {
    min = atof(     argv[9] );
    max = atof(     argv[10] );
    }
	

	// Graphics Factory
	vtkGraphicsFactory * graphics_factory 
		= vtkGraphicsFactory::New();
	graphics_factory->SetOffScreenOnlyMode( 1);
	graphics_factory->SetUseMesaClasses( 1 );

	// Imaging Factory
	vtkImagingFactory * imaging_factory 
		= vtkImagingFactory::New();
	imaging_factory->SetUseMesaClasses( 1 );

	// Image Reader
	vtkMetaImageReader * reader 
		= vtkMetaImageReader::New();
	reader->SetFileName( in_file_name );

	// Get center
	reader->Update();
	vtkImageData * image = reader->GetOutput();
  double * range = image->GetScalarRange();
  if ( argc > 9 )
    {
    range[0] = min;
    range[1] = max;
    }
	double * image_spacing = image->GetSpacing();
	double * image_origin = image->GetOrigin();
	int * whole_extent = image->GetWholeExtent();
	int image_dim[3];
	image_dim[0] = whole_extent[1] + 1;
	image_dim[1] = whole_extent[3] + 1;
	image_dim[2] = whole_extent[5] + 1;
	double fpx = image_origin[0] 
		+ image_dim[0] / 2 * image_spacing[0];
	double fpy = image_origin[1] 
		+ image_dim[1] / 2 * image_spacing[1];
	double fpz = image_origin[2] 
		+ image_dim[2] / 2 * image_spacing[2];

  // Pixel Diagonal
  double pix_diag = sqrt(
      image_spacing[0] * image_spacing[0] +
      image_spacing[1] * image_spacing[1] +
      image_spacing[2] * image_spacing[2] );
	
	// Opacity
	vtkPiecewiseFunction * opacity 
		= vtkPiecewiseFunction::New();
  opacity->AddPoint( range[0], 0.0 );
  opacity->AddPoint( range[1], 1.0 );
    

	// Color
	vtkColorTransferFunction * color 
		= vtkColorTransferFunction::New();
  color->AddRGBPoint( range[0],                                       0.0, 0.0, 0.0 );
  color->AddRGBPoint( range[0] + 1.0 * ( range[1] - range[0] ) / 3.0, 1.0, 0.0, 0.0 );
  color->AddRGBPoint( range[0] + 2.0 * ( range[1] - range[0] ) / 3.0, 1.0, 1.0, 0.0 );
  color->AddRGBPoint( range[1],                                       1.0, 1.0, 1.0 );
    
 
	// Volume Property
	vtkVolumeProperty * vol_prop = vtkVolumeProperty::New();
	vol_prop->SetColor( color );
	vol_prop->SetScalarOpacity( opacity );
	//vol_prop->ShadeOff();
	//vol_prop->SetInterpolationTypeToLinear();

	// Mapper
	vtkFixedPointVolumeRayCastMapper * mapper
		= vtkFixedPointVolumeRayCastMapper::New();
  mapper->SetImageSampleDistance( 0.25 );
  mapper->SetSampleDistance( pix_diag / 5.0 );
	mapper->SetInputConnection( reader->GetOutputPort() );

	// Volume
	vtkVolume * volume = vtkVolume::New();
	volume->SetMapper( mapper );
	volume->SetProperty( vol_prop );

	// Camera
	vtkCamera * camera = vtkCamera::New();
	camera->SetViewUp (vux, vuy, vuz);
	camera->SetPosition (
			fpx + px, 
			fpy + py, 
			fpz + pz );
	camera->SetFocalPoint (fpx, fpy, fpz);
	camera->ComputeViewPlaneNormal();

	// Renderer
	vtkRenderer * renderer = vtkRenderer::New();
	renderer->SetBackground( 0.4, 0.4, 0.4 );
  renderer->AddVolume( volume );
	renderer->SetActiveCamera(camera);
	renderer->ResetCamera();

	//camera->Zoom(3.0);
	
	// Render Window
	vtkRenderWindow * render_window = vtkRenderWindow::New();
	render_window->SetSize( WSIZE, WSIZE );
	render_window->SetOffScreenRendering( 1 );
	render_window->AddRenderer( renderer );
	// ***** STEREO *****
	//render_window->StereoRenderOn();

	// Window to Image	
	vtkWindowToImageFilter * win_2_image 
		= vtkWindowToImageFilter::New();
	win_2_image->SetInput( render_window );

	// PNG Writer
	vtkPNGWriter * writer = vtkPNGWriter::New();
	writer->SetFileName( out_file_name );
	writer->SetInputConnection( win_2_image->GetOutputPort() );

	render_window->Render();

	writer->Write();
	
	return 0;
}

