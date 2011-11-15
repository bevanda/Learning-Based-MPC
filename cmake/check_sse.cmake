# borrowed from:
# https://code.ros.org/trac/ros/ticket/2594

# check for SSE extensions
include(CheckCXXSourceRuns)
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
 set(SSE_FLAGS)

 set(CMAKE_REQUIRED_FLAGS "-msse3")
 check_cxx_source_runs("
  #include <pmmintrin.h>

  int main()
  {
     __m128d a, b;
     double vals[2] = {0};
     a = _mm_loadu_pd(vals);
     b = _mm_hadd_pd(a,a);
     _mm_storeu_pd(vals, b);
     return 0;
  }"
  HAS_SSE3_EXTENSIONS)

 set(CMAKE_REQUIRED_FLAGS "-msse2")
 check_cxx_source_runs("
  #include <emmintrin.h>

  int main()
  {
      __m128d a, b;
      double vals[2] = {0};
      a = _mm_loadu_pd(vals);
      b = _mm_add_pd(a,a);
      _mm_storeu_pd(vals,b);
      return 0;
   }"
   HAS_SSE2_EXTENSIONS)

 set(CMAKE_REQUIRED_FLAGS "-msse")
 check_cxx_source_runs("
  #include <xmmintrin.h>
  int main()
  {
      __m128 a, b;
      float vals[4] = {0};
      a = _mm_loadu_ps(vals);
      b = a;
      b = _mm_add_ps(a,b);
      _mm_storeu_ps(vals,b);
      return 0;
  }"
  HAS_SSE_EXTENSIONS)

 set(CMAKE_REQUIRED_FLAGS)

 if(HAS_SSE3_EXTENSIONS)
  message(STATUS "Using SSE3 extensions")
  set(SSE_FLAGS "-msse3 -mfpmath=sse")
 elseif(HAS_SSE2_EXTENSIONS)
  message(STATUS "Using SSE2 extensions")
  set(SSE_FLAGS "-msse2 -mfpmath=sse")
 elseif(HAS_SSE_EXTENSIONS)
  message(STATUS "Using SSE extensions")
  set(SSE_FLAGS "-msse -mfpmath=sse")
 endif()

 add_definitions(${SSE_FLAGS})
elseif(MSVC)
 check_cxx_source_runs("
  #include <emmintrin.h>

  int main()
  {
      __m128d a, b;
      double vals[2] = {0};
      a = _mm_loadu_pd(vals);
      b = _mm_add_pd(a,a);
      _mm_storeu_pd(vals,b);
      return 0;
   }"
   HAS_SSE2_EXTENSIONS)
 if( HAS_SSE2_EXTENSIONS )
  message(STATUS "Using SSE2 extensions")
  add_definitions( "/arch:SSE2 /fp:fast -D__SSE__ -D__SSE2__" )
 endif()
endif()
