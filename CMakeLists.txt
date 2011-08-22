INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#Change PROJECT_NAME to the name of your project
PROJECT(MPIFastMarching)

cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
ITK
VTK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

INCLUDE( FindMPI )
INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(ITK_Libraries
ITKCommon
ITKIO
)

SET(MPI_Libraries
"${MPI_LIBRARY}"
"${MPI_EXTRA_LIBRARY}"
)

SET(VTK_Libraries
vtkImaging
vtkIO
vtkCommon
vtkRendering
vtkVolumeRendering
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "MPIFastMarching")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${ITK_Libraries} ${MPI_Libraries})

SET(CurrentExe "Difference")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${ITK_Libraries})

SET(CurrentExe "OSVolumeRender")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${VTK_Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
ADD_TEST(TestVolumeRenderInput ./OSVolumeRender
  ${CMAKE_SOURCE_DIR}/input.mhd input_top_test.png 0 1 0 0 0 1)


ADD_TEST(TestSerialFastMarching mpirun -n 1
  ./MPIFastMarching ${CMAKE_SOURCE_DIR}/input.mhd
  ${CMAKE_SOURCE_DIR}/seeds.txt fast_n1_i0_test.mhd 0)


ADD_TEST(TestVolumeRenderSerialFastMarching
  ./OSVolumeRender fast_n1_i0_test.mhd fast_n1_i0_top_test.png
  0 1 0  0 0  1  200 0)


ADD_TEST(TestMPIFastMarching mpirun -n 4
  ./MPIFastMarching
  ${CMAKE_SOURCE_DIR}/input.mhd ${CMAKE_SOURCE_DIR}/seeds.txt
  fast_n4_i8_test.mhd 8)


ADD_TEST(TestVolumeRenderMPIFastMarching
  ./OSVolumeRender fast_n4_i8_test.mhd fast_n4_i8_top_test.png
  0 1 0  0 0  1  200 0)


ADD_TEST(TestDifference ./Difference fast_n4_i8_test.mhd
  fast_n1_i0_test.mhd diff_test.mhd)


ADD_TEST(TestVolumeRenderDifference ./OSVolumeRender
  diff_test.mhd diff_top_test.png  0 1 0  0 0 1)


ADD_TEST(TestCompareVolumeRenderInput ImageCompare input_top_test.png
  ${CMAKE_SOURCE_DIR}/input_top.png)


ADD_TEST(TestCompareSerialFastMarching ImageCompare fast_n1_i0_test.mhd
  ${CMAKE_SOURCE_DIR}/fast_n1_i0.mhd)


ADD_TEST(TestCompareVolumeRenderSerialFastMarching ImageCompare
  fast_n1_i0_top_test.png
  ${CMAKE_SOURCE_DIR}/fast_n1_i0_top.png)


ADD_TEST(TestCompareMPIFastMarching ImageCompare fast_n4_i8_test.mhd
  ${CMAKE_SOURCE_DIR}/fast_n4_i8.mhd)


ADD_TEST(TestCompareVolumeRenderMPIFastMarching ImageCompare
  fast_n4_i8_top_test.png
  ${CMAKE_SOURCE_DIR}/fast_n4_i8_top.png)


ADD_TEST(TestCompareDifference ImageCompare diff_test.mhd ${CMAKE_SOURCE_DIR}/diff.mhd)


ADD_TEST(TestCompareVolumeRenderDifference ImageCompare diff_top_test.png
  ${CMAKE_SOURCE_DIR}/diff_top.png)

# This test fails.
ADD_TEST(TestCompareSerialMPI ImageCompare fast_n4_i8_test.mhd fast_n1_i0_test.mhd)

#Once you're done modifying this template, you should rename it to "CMakeLists.txt"
