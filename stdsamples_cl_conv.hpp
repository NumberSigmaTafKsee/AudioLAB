#include <complex>
#include <iostream>
#include <string>

#ifdef __MACH__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace cl_conv {

  inline const char *cl_string(int err) {
    switch (err) {
    case CL_SUCCESS:
      return "Success!";
    case CL_DEVICE_NOT_FOUND:
      return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE:
      return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE:
      return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
      return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES:
      return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY:
      return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE:
      return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP:
      return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH:
      return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:
      return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE:
      return "Program build failure";
    case CL_MAP_FAILURE:
      return "Map failure";
    case CL_INVALID_VALUE:
      return "Invalid value";
    case CL_INVALID_DEVICE_TYPE:
      return "Invalid device type";
    case CL_INVALID_PLATFORM:
      return "Invalid platform";
    case CL_INVALID_DEVICE:
      return "Invalid device";
    case CL_INVALID_CONTEXT:
      return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES:
      return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE:
      return "Invalid command queue";
    case CL_INVALID_HOST_PTR:
      return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT:
      return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
      return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE:
      return "Invalid image size";
    case CL_INVALID_SAMPLER:
      return "Invalid sampler";
    case CL_INVALID_BINARY:
      return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS:
      return "Invalid build options";
    case CL_INVALID_PROGRAM:
      return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE:
      return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME:
      return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION:
      return "Invalid kernel definition";
    case CL_INVALID_KERNEL:
      return "Invalid kernel";
    case CL_INVALID_ARG_INDEX:
      return "Invalid argument index";
    case CL_INVALID_ARG_VALUE:
      return "Invalid argument value";
    case CL_INVALID_ARG_SIZE:
      return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS:
      return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION:
      return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE:
      return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE:
      return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET:
      return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST:
      return "Invalid event wait list";
    case CL_INVALID_EVENT:
      return "Invalid event";
    case CL_INVALID_OPERATION:
      return "Invalid operation";
    case CL_INVALID_GL_OBJECT:
      return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE:
      return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL:
      return "Invalid mip-map level";
    default:
      return "Unknown error";
    }
  }

  class Clpconv {
    int N, bins;
    int bsize, nparts, wp, wp2;
    cl_mem w[2], w2[2], b;
    cl_mem in1, in2, out;
    cl_mem spec1, spec2, olap;
    cl_context context;
    cl_command_queue commands1, commands2;
    cl_program program;
    cl_kernel fft_kernel1, reorder_kernel1;
    cl_kernel fft_kernel2, reorder_kernel2;
    cl_kernel r2c_kernel1, r2c_kernel2, c2r_kernel;
    cl_kernel convol_kernel, olap_kernel;
    void (*err)(std::string s, void *uData);
    void *userData;
    int cl_err;
    int host_mem;

    static void msg(std::string str, void *userData) {
      if (userData == NULL)
        std::cout << str << std::endl;
    }

  public:
    /** Constructor \n
        device_id - OpenCL device ID \n
        cvs - impulse response size \n
        pts - partition size \n
        errs - error message callback \n
        uData - callback user data \n
    */
    Clpconv(cl_device_id device_id, int cvs, int pts,
            void (*errs)(std::string s, void *d) = NULL, void *uData = NULL,
            void *in1 = NULL, void *in2 = NULL, void *out = 0);
    ~Clpconv();

    /** returns an error string relative to error code err */
    const char *cl_error_string(int err) { return cl_string(err); }

    /** set the convolution impulse response
        ir - impulse response of size cvs;
    */
    int push_ir(float *ir);

    /** Convolution computation
        output - output array (partition size samples) \n
        input - input array (partition size * 2, but
          holding only partition size samples, zero-padded to
          partition length)  \n
    */
    int convolution(float *output, float *input);

    /** Time-varying convolution computation
        output - output array (partition size samples) \n
        input1, input2 - input arrays (partition size * 2, but
          holding only partition size samples, zero-padded to
          partition length)  \n
    */
    int convolution(float *output, float *input1, float *input2);

    /** get a recorded error code, CL_SUCCESS if no error was recorded
     */
    int get_cl_err() { return cl_err; }
  };

  class Cldconv {
    int irsize, vsize, wp;
    cl_mem buff, coefs, del;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel convol_kernel;
    size_t wgs;
    void (*err)(std::string s, void *uData);
    void *userData;
    int cl_err;

    static void msg(std::string str, void *userData) {
        if (userData == NULL)
        std::cout << str << std::endl;
    }

    public:
    /** Constructor \n
         device_id - OpenCL device ID \n
        cvs - impulse response size \n
        vsize - processing vector size \n
        errs - error message callback \n
        uData - callback user data \n
    */
    Cldconv(cl_device_id device_id, int cvs, int vsize,
            void (*errs)(std::string s, void *d) = NULL, void *uData = NULL);

    ~Cldconv();

    /** returns an error string relative to error code err */
    const char *cl_error_string(int err) { return cl_string(err); }

    /** set the convolution impulse response
         ir - impulse response of size cvs;
    */
    int push_ir(float *ir);

    /** Convolution computation
         output - output array (vsize samples) \n
        input - input array (vsize samples) \n
    */
    int convolution(float *output, float *input);

    int convolution(float *out, float *in1, float *in2);

    /** get a recorded error code, CL_SUCCESS if no error was recorded
     */
    int get_cl_err() { return cl_err; }
    };
  const double PI = 3.141592653589793;

  const char *pconvcode = R"(
  /* complex type */
  typedef float2 cmplx;
  /* complex product */
  inline cmplx prod(cmplx a, cmplx b){
      return (cmplx)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); 
  }
  /* complex conj */
  inline cmplx conjg(cmplx a) {
      return (cmplx) (a.x, - a.y);
  }
  /* rotation by pi */
  inline cmplx rot(cmplx a) {
    return (cmplx) (-a.y, a.x);
  } 
  /* atomic add */
  inline void AtomicAdd(volatile __global float *source, const float operand) {
      union {
          uint intVal;
          float floatVal;
      } newVal;
      union {
          uint intVal;
          float floatVal;
      } prevVal;
      do {
          prevVal.floatVal = *source;
          newVal.floatVal = prevVal.floatVal + operand;
      } while (atomic_cmpxchg((volatile __global uint *) source, 
        prevVal.intVal, newVal.intVal) != prevVal.intVal);
  }
  /* data reordering */
  kernel void reorder(global cmplx *out, global cmplx *in, global const int *b, int offs) {
    int k = get_global_id(0);
    out += offs;
    out[k] = in[b[k]]; 
    in[b[k]] = 0.f;
  }
  /* fft stage  */
  kernel void fft(global cmplx *s, global const cmplx *w, int N, int n2, int offs) {
  int k, i, m, n;
  cmplx e, o;
  s += offs;
  k = get_global_id(0)*n2;
  m = k/N; 
  n = n2 >> 1;
  k =  k%N + m;
  i = k + n;
  e = s[k];
  o = prod(s[i],w[m*N/n2]);
  s[k] = e + o;
  s[i] = e - o;  
  }
  /* rfft conversion */
  kernel void r2c(global cmplx *c, global const cmplx *w, int N, int offs) {
    int i = get_global_id(0);
    if(!i) {
    c[0] = (cmplx) ((c[0].x + c[0].y)*.5f, (c[0].x - c[0].y)*.5f);
    return;
    }
    int j = N - i;
    cmplx e, o, cj = conjg(c[j]), p;
    c += offs;
    e = .5f*(c[i] + cj);
    o = .5f*rot(cj - c[i]);
    p = prod(w[i], o); 
    c[i] = e + p;
    c[j] = conjg(e - p);
  }
  /* inverse rfft conversion */
  kernel void c2r(global cmplx *c, global const cmplx *w, int N) {
    int i = get_global_id(0);
    if(!i) {
    c[0] = (cmplx) ((c[0].x + c[0].y), (c[0].x - c[0].y));
    return; 
    }
    int j = N - i;
    cmplx e, o, cj = conjg(c[j]), p;
    e = .5f*(c[i] + cj);
    o = .5f*rot(c[i] - cj);
    p = prod(w[i], o);
    c[i] = e + p;
    c[j] = conjg(e - p); 
  }
  /* convolution */
  kernel void convol(global float *out, global const cmplx *in, 
          global const cmplx *coef, int rp, int b, 
          int nparts) {
    /* thread count */
    int k = get_global_id(0); /* bin pos */ 
    int n = k%b;  /* inframe pos   */
    int n2 = n << 1;
    cmplx s;
    rp += k/b;       /*  rp pos */
    /* select correct input buffer */
    in += (rp < nparts ? rp : rp - nparts)*b;  
    /* complex multiplication + sums */
    s = n ? prod(in[n], coef[k]) : 
              (cmplx) (in[0].x*coef[k].x, in[0].y*coef[k].y);
    AtomicAdd(&out[n2], s.x);  
    AtomicAdd(&out[n2+1], s.y);                                
  }  
  /* sample-by-sample overlap-add operation */
  kernel void olap(global float *buf, global const float *in, int parts){
    int n = get_global_id(0);
    buf[n] = (in[n] + buf[parts+n])/parts;
    buf[parts+n] = in[parts+n];
  }
  )";

  /*  kernel dispatch functions */
  inline int reorder(cl_mem *out, cl_mem *in, cl_mem *b, int offs,
                    cl_command_queue commands, cl_kernel kern, size_t threads) {
    clSetKernelArg(kern, 3, sizeof(cl_int), &offs);
    clSetKernelArg(kern, 2, sizeof(cl_mem), b);
    clSetKernelArg(kern, 1, sizeof(cl_mem), in);
    clSetKernelArg(kern, 0, sizeof(cl_mem), out);
    return clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                  NULL, NULL);
  }

  inline int fft(cl_mem *data, cl_mem *w, int bins, int offs,
                cl_command_queue commands, cl_kernel kern, size_t threads) {
    int cl_err, n2;
    for (int n = 1; n < bins; n *= 2) {
      n2 = n << 1;
      clSetKernelArg(kern, 4, sizeof(cl_int), &offs);
      clSetKernelArg(kern, 3, sizeof(cl_int), &n2);
      clSetKernelArg(kern, 2, sizeof(cl_int), &bins);
      clSetKernelArg(kern, 1, sizeof(cl_mem), w);
      clSetKernelArg(kern, 0, sizeof(cl_mem), data);
      cl_err = clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                      NULL, NULL);
    }
    return cl_err;
  }

  inline int real_cmplx(cl_mem *data, cl_mem *w, int bins, int offs,
                        cl_command_queue commands, cl_kernel kern,
                        size_t threads) {
    clSetKernelArg(kern, 3, sizeof(cl_int), &offs);
    clSetKernelArg(kern, 2, sizeof(cl_int), &bins);
    clSetKernelArg(kern, 1, sizeof(cl_mem), w);
    clSetKernelArg(kern, 0, sizeof(cl_mem), data);
    return clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                  NULL, NULL);
  }

  inline int cmplx_real(cl_mem *data, cl_mem *w, int bins,
                        cl_command_queue commands, cl_kernel kern,
                        size_t threads) {
    clSetKernelArg(kern, 2, sizeof(cl_int), &bins);
    clSetKernelArg(kern, 1, sizeof(cl_mem), w);
    clSetKernelArg(kern, 0, sizeof(cl_mem), data);
    return clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                  NULL, NULL);
  }

  inline int convol(cl_mem *out, cl_mem *in, cl_mem *coefs, int wp, int bins,
                    int nparts, cl_command_queue commands, cl_kernel kern,
                    size_t threads) {
    clSetKernelArg(kern, 5, sizeof(cl_int), &nparts);
    clSetKernelArg(kern, 4, sizeof(cl_int), &bins);
    clSetKernelArg(kern, 3, sizeof(cl_int), &wp);
    clSetKernelArg(kern, 2, sizeof(cl_mem), coefs);
    clSetKernelArg(kern, 1, sizeof(cl_mem), in);
    clSetKernelArg(kern, 0, sizeof(cl_mem), out);
    return clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                  NULL, NULL);
  }

  inline int ola(cl_mem *out, cl_mem *in, int parts, cl_command_queue commands,
                cl_kernel kern, size_t threads) {
    clSetKernelArg(kern, 2, sizeof(cl_int), &parts);
    clSetKernelArg(kern, 1, sizeof(cl_mem), in);
    clSetKernelArg(kern, 0, sizeof(cl_mem), out);
    return clEnqueueNDRangeKernel(commands, kern, 1, NULL, &threads, NULL, 0,
                                  NULL, NULL);
  }

  Clpconv::Clpconv(cl_device_id device_id, int cvs, int pts,
                  void (*errs)(std::string s, void *d), void *uData, void *inp1,
                  void *inp2, void *outp)
      : N(pts << 1), bins(pts), bsize((cvs / pts) * bins), nparts(cvs / pts),
        wp(0), wp2(nparts - 1), w{NULL, NULL}, w2{NULL, NULL}, b(NULL), in1(NULL),
        in2(NULL), out(NULL), spec1(NULL), spec2(NULL), olap(NULL), context(NULL),
        commands1(NULL), commands2(NULL), program(NULL), fft_kernel1(NULL),
        reorder_kernel1(NULL), fft_kernel2(NULL), reorder_kernel2(NULL),
        r2c_kernel1(NULL), r2c_kernel2(NULL), c2r_kernel(NULL),
        convol_kernel(NULL), olap_kernel(NULL),
        err(errs == NULL ? this->msg : errs), userData(uData), cl_err(CL_SUCCESS),
        host_mem(((uintptr_t)inp1 & (uintptr_t)inp2 & (uintptr_t)out) ? 1 : 0) {

    context = clCreateContext(0, 1, &device_id, NULL, NULL, &cl_err);
    if (!context) {
      err(cl_error_string(cl_err), userData);
      return;
    }

    /* two command queues for task parallelism */
    commands1 = clCreateCommandQueue(context, device_id, 0, &cl_err);
    if (!commands1) {
      err(cl_error_string(cl_err), userData);
      return;
    }
    commands2 = clCreateCommandQueue(context, device_id, 0, &cl_err);
    if (!commands2) {
      err(cl_error_string(cl_err), userData);
      return;
    }

    program = clCreateProgramWithSource(context, 1, (const char **)&pconvcode,
                                        NULL, &cl_err);
    if (!program) {
      err("error creating opencl program\n", userData);
      err(cl_error_string(cl_err), userData);
      return;
    }
    cl_err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (cl_err) {
      char log[2048];
      size_t llen;
      clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(log),
                            log, &llen);
      err("error building opencl program\n", userData);
      err(cl_error_string(cl_err), userData);
      err(log, userData);
      return;
    }

    /* separate kernel objects for each command queue*/
    reorder_kernel1 = clCreateKernel(program, "reorder", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    fft_kernel1 = clCreateKernel(program, "fft", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    r2c_kernel1 = clCreateKernel(program, "r2c", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    reorder_kernel2 = clCreateKernel(program, "reorder", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    fft_kernel2 = clCreateKernel(program, "fft", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    r2c_kernel2 = clCreateKernel(program, "r2c", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    c2r_kernel = clCreateKernel(program, "c2r", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    convol_kernel = clCreateKernel(program, "convol", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);
    olap_kernel = clCreateKernel(program, "olap", &cl_err);
    if (cl_err != 0)
      err(cl_error_string(cl_err), userData);

    in1 = clCreateBuffer(context, in1 ? CL_MEM_USE_HOST_PTR : 0,
                        bins * sizeof(cl_float2), inp1, &cl_err);
    in2 = clCreateBuffer(context, in2 ? CL_MEM_USE_HOST_PTR : 0,
                        bins * sizeof(cl_float2), inp2, &cl_err);
    olap = clCreateBuffer(context, out ? CL_MEM_USE_HOST_PTR : 0,
                          bins * sizeof(cl_float2), outp, &cl_err);
    out = clCreateBuffer(context, 0, bins * sizeof(cl_float2), NULL, &cl_err);
    spec1 = clCreateBuffer(context, CL_MEM_READ_ONLY, bsize * sizeof(cl_float2),
                          NULL, &cl_err);
    spec2 = clCreateBuffer(context, CL_MEM_READ_ONLY, bsize * sizeof(cl_float2),
                          NULL, &cl_err);
    w[0] = clCreateBuffer(context, CL_MEM_READ_ONLY, bins * sizeof(cl_float2),
                          NULL, &cl_err);
    w[1] = clCreateBuffer(context, CL_MEM_READ_ONLY, bins * sizeof(cl_float2),
                          NULL, &cl_err);
    w2[0] = clCreateBuffer(context, CL_MEM_READ_ONLY, bins * sizeof(cl_float2),
                          NULL, &cl_err);
    w2[1] = clCreateBuffer(context, CL_MEM_READ_ONLY, bins * sizeof(cl_float2),
                          NULL, &cl_err);
    b = clCreateBuffer(context, CL_MEM_READ_ONLY, bins * sizeof(cl_int), NULL,
                      NULL);
    /* twiddle */
    std::vector<std::complex<float>> wd(bins);
    for (int i = 0; i < bins; i++) {
      wd[i].real(cos(i * 2 * PI / bins));
      wd[i].imag(-sin(i * 2 * PI / bins));
    }
    clEnqueueWriteBuffer(commands1, w[0], CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)wd.data(), 0, NULL, NULL);
    for (int i = 0; i < bins; i++) {
      wd[i].real(cos(i * 2 * PI / bins));
      wd[i].imag(sin(i * 2 * PI / bins));
    }
    clEnqueueWriteBuffer(commands1, w[1], CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)wd.data(), 0, NULL, NULL);

    for (int i = 0; i < bins; i++) {
      wd[i].real(cos(i * PI / bins));
      wd[i].imag(-sin(i * PI / bins));
    }
    clEnqueueWriteBuffer(commands1, w2[0], CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)wd.data(), 0, NULL, NULL);
    for (int i = 0; i < bins; i++) {
      wd[i].real(cos(i * PI / bins));
      wd[i].imag(sin(i * PI / bins));
    }
    clEnqueueWriteBuffer(commands1, w2[1], CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)wd.data(), 0, NULL, NULL);

    /* bit-reversed indices */
    std::vector<int> bp(bins);
    for (int i = 0; i < bins; i++)
      bp[i] = i;
    for (int i = 1, n = bins / 2; i<bins; i = i << 1, n = n>> 1)
      for (int j = 0; j < i; j++)
        bp[i + j] = bp[j] + n;
    clEnqueueWriteBuffer(commands1, b, CL_TRUE, 0, sizeof(cl_int) * bins,
                        (const void *)bp.data(), 0, NULL, NULL);

    std::vector<float> zeros(bsize * 2, 0);
    clEnqueueWriteBuffer(commands1, olap, CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)zeros.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(commands1, spec1, CL_TRUE, 0, sizeof(cl_float2) * bsize,
                        (const void *)zeros.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(commands1, spec2, CL_TRUE, 0, sizeof(cl_float2) * bsize,
                        (const void *)zeros.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(commands1, in1, CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)zeros.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(commands1, in2, CL_TRUE, 0, sizeof(cl_float2) * bins,
                        (const void *)zeros.data(), 0, NULL, NULL);
  }

  Clpconv::~Clpconv() {
    clReleaseMemObject(w2[0]);
    clReleaseMemObject(w2[1]);
    clReleaseMemObject(w[0]);
    clReleaseMemObject(w[1]);
    clReleaseMemObject(b);
    clReleaseMemObject(out);
    clReleaseMemObject(in2);
    clReleaseMemObject(in1);
    clReleaseMemObject(olap);
    clReleaseMemObject(spec1);
    clReleaseMemObject(spec2);
    clReleaseKernel(fft_kernel2);
    clReleaseKernel(reorder_kernel2);
    clReleaseKernel(r2c_kernel2);
    clReleaseKernel(fft_kernel1);
    clReleaseKernel(reorder_kernel1);
    clReleaseKernel(r2c_kernel1);
    clReleaseKernel(c2r_kernel);
    clReleaseKernel(olap_kernel);
    clReleaseKernel(convol_kernel);
    clReleaseCommandQueue(commands1);
    clReleaseCommandQueue(commands2);
    clReleaseContext(context);
  }

  int Clpconv::push_ir(float *ir) {
    size_t bytes = sizeof(cl_float2) * bins;
    size_t threads;
    int n2;
    for (int i = 0; i < nparts; i++) {
      if (!host_mem)
        cl_err = clEnqueueWriteBuffer(commands1, in2, CL_TRUE, 0, bytes >> 1,
                                      &ir[i * bins], 0, NULL, NULL);
      if (cl_err != CL_SUCCESS)
        return cl_err;
      cl_err =
          reorder(&spec2, &in2, &b, wp2 * bins, commands1, reorder_kernel1, bins);
      if (cl_err != CL_SUCCESS)
        return cl_err;
      cl_err =
          fft(&spec2, &w[0], bins, wp2 * bins, commands1, fft_kernel1, bins >> 1);
      if (cl_err != CL_SUCCESS)
        return cl_err;
      cl_err = real_cmplx(&spec2, &w2[0], bins, wp2 * bins, commands1,
                          r2c_kernel1, bins >> 1);
      if (cl_err != CL_SUCCESS)
        return cl_err;
      wp2 = wp2 == 0 ? nparts - 1 : wp2 - 1;
    }
    return CL_SUCCESS;
  }

  int Clpconv::convolution(float *output, float *input) {
    int n2;
    size_t bytes = sizeof(cl_float2) * bins;
    char zro = 0;
    if (!host_mem)
      cl_err = clEnqueueWriteBuffer(commands1, in1, 0, 0, bytes >> 1, input, 0,
                                    NULL, NULL);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        reorder(&spec1, &in1, &b, wp * bins, commands1, reorder_kernel1, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        fft(&spec1, &w[0], bins, wp * bins, commands1, fft_kernel1, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = real_cmplx(&spec1, &w2[0], bins, wp * bins, commands1, r2c_kernel1,
                        bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    wp = wp != nparts - 1 ? wp + 1 : 0;
    cl_err = convol(&in1, &spec1, &spec2, wp, bins, nparts, commands1,
                    convol_kernel, bsize);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = cmplx_real(&in1, &w2[1], bins, commands1, c2r_kernel, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = reorder(&out, &in1, &b, 0, commands1, reorder_kernel1, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = fft(&out, &w[1], bins, 0, commands1, fft_kernel1, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = ola(&olap, &out, bins, commands1, olap_kernel, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    if (!host_mem)
      cl_err = clEnqueueReadBuffer(commands1, olap, CL_TRUE, 0, bytes >> 1,
                                  output, 0, NULL, NULL);
    return cl_err;
  }

  int Clpconv::convolution(float *output, float *input1, float *input2) {
    size_t bytes = sizeof(cl_float2) * bins;
    if (!host_mem)
      cl_err = clEnqueueWriteBuffer(commands1, in1, 0, 0, bytes >> 1, input1, 0,
                                    NULL, NULL);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    if (!host_mem)
      cl_err = clEnqueueWriteBuffer(commands2, in2, 0, 0, bytes >> 1, input2, 0,
                                    NULL, NULL);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        reorder(&spec1, &in1, &b, wp * bins, commands1, reorder_kernel1, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        reorder(&spec2, &in2, &b, wp2 * bins, commands2, reorder_kernel2, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        fft(&spec1, &w[0], bins, wp * bins, commands1, fft_kernel1, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err =
        fft(&spec2, &w[0], bins, wp2 * bins, commands2, fft_kernel2, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = real_cmplx(&spec1, &w2[0], bins, wp * bins, commands1, r2c_kernel1,
                        bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = real_cmplx(&spec2, &w2[0], bins, wp2 * bins, commands2, r2c_kernel2,
                        bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    wp = wp != nparts - 1 ? wp + 1 : 0;
    wp2 = wp2 == 0 ? nparts - 1 : wp2 - 1;
    clFinish(commands2);
    cl_err = convol(&in1, &spec1, &spec2, wp, bins, nparts, commands1,
                    convol_kernel, bsize);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = cmplx_real(&in1, &w2[1], bins, commands1, c2r_kernel, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = reorder(&out, &in1, &b, 0, commands1, reorder_kernel1, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = fft(&out, &w[1], bins, 0, commands1, fft_kernel1, bins >> 1);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    cl_err = ola(&olap, &out, bins, commands1, olap_kernel, bins);
    if (cl_err != CL_SUCCESS)
      return cl_err;
    if (!host_mem)
      cl_err = clEnqueueReadBuffer(commands1, olap, CL_TRUE, 0, bytes >> 1,
                                  output, 0, NULL, NULL);
    return cl_err;
  }

  const char *dconvcode = R"(
    /* atomic add */
    inline void AtomicAdd(volatile __global float *source, const float operand) {
        union {
            uint intVal;
            float floatVal;
        } newVal;
        union {
            uint intVal;
            float floatVal;
        } prevVal;
        do {
            prevVal.floatVal = *source;
            newVal.floatVal = prevVal.floatVal + operand;
        } while (atomic_cmpxchg((volatile __global uint *) source, 
        prevVal.intVal, newVal.intVal) != prevVal.intVal);
    }
    kernel void convol(global float *out, global const float *del, global const 
            float *coefs, int irsize, int rp, int vsize) {
    int t = get_global_id(0);
    float tap;
    if(t >= irsize*vsize) return;
    int n =  t%vsize;  /* sample index */
    int h =  t/vsize;  /* coeff index */
    int end = irsize+vsize;
    rp += n + h; /* read point, oldest -> newest */
    tap = del[rp < end ? rp : rp%end]*coefs[irsize-1-h];  /* single tap */
    AtomicAdd(&out[n], tap);
    }
    )";

    Cldconv::Cldconv(cl_device_id device_id, int cvs, int vsiz,
                    void (*errs)(std::string s, void *d), void *uData)
        : irsize(cvs), vsize(vsiz), wp(0), buff(NULL), coefs(NULL), del(NULL),
        context(NULL), commands(NULL), program(NULL), convol_kernel(NULL), wgs(0),
        err(errs == NULL ? this->msg : errs), userData(uData),
        cl_err(CL_SUCCESS) {

    context = clCreateContext(0, 1, &device_id, NULL, NULL, &cl_err);
    if (!context) {
        err(cl_error_string(cl_err), userData);
        return;
    }
    commands = clCreateCommandQueue(context, device_id, 0, &cl_err);
    if (!commands) {
        err(cl_error_string(cl_err), userData);
        return;
    }
    program = clCreateProgramWithSource(context, 1, (const char **)&dconvcode,
                                        NULL, &cl_err);
    if (!program) {
        err("error creating conv program\n", userData);
        err(cl_error_string(cl_err), userData);
        return;
    }
    cl_err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (cl_err) {
        char log[2048];
        size_t llen;
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(log),
                            log, &llen);
        err("error building conv program\n", userData);
        err(cl_error_string(cl_err), userData);
        err(log, userData);
        return;
    }
    convol_kernel = clCreateKernel(program, "convol", &cl_err);
    clGetKernelWorkGroupInfo(convol_kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE,
                            sizeof(wgs), &wgs, NULL);
    if (wgs > vsize * irsize)
        wgs = vsize * irsize;

    buff = clCreateBuffer(context, 0, vsize * sizeof(cl_float), NULL, NULL);
    del = clCreateBuffer(context, CL_MEM_READ_ONLY,
                        (irsize + vsize) * sizeof(cl_float), NULL, NULL);
    coefs = clCreateBuffer(context, CL_MEM_READ_ONLY,
                            (irsize + vsize) * sizeof(cl_float), NULL, NULL);

    clSetKernelArg(convol_kernel, 0, sizeof(cl_mem), &buff);
    clSetKernelArg(convol_kernel, 1, sizeof(cl_mem), &del);
    clSetKernelArg(convol_kernel, 2, sizeof(cl_mem), &coefs);
    clSetKernelArg(convol_kernel, 3, sizeof(cl_int), &irsize);
    clSetKernelArg(convol_kernel, 5, sizeof(cl_int), &vsize);
    }

    Cldconv::~Cldconv() {
    clReleaseMemObject(del);
    clReleaseMemObject(buff);
    clReleaseMemObject(coefs);
    clReleaseKernel(convol_kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    }

    int Cldconv::convolution(float *out, float *in) {
    size_t bytes = vsize * sizeof(cl_float), threads = irsize * vsize;
    char zro = 0;
    if (wp > irsize) {
        int front = wp - irsize;
        bytes = (vsize - front) * sizeof(cl_float);
        clEnqueueWriteBuffer(commands, del, CL_TRUE, wp * sizeof(cl_float), bytes,
                            in, 0, NULL, NULL);
        bytes = front * sizeof(cl_float);
        clEnqueueWriteBuffer(commands, del, CL_TRUE, 0, bytes, &in[vsize - front],
                            0, NULL, NULL);
    } else
        clEnqueueWriteBuffer(commands, del, CL_TRUE, wp * sizeof(cl_float), bytes,
                            in, 0, NULL, NULL);
    clEnqueueFillBuffer(commands, buff, &zro, 1, 0, bytes, 0, NULL, NULL);
    wp = (wp + vsize) % (irsize + vsize);
    clSetKernelArg(convol_kernel, 4, sizeof(cl_int), &wp);
    cl_err = clEnqueueNDRangeKernel(commands, convol_kernel, 1, NULL, &threads,
                                    &wgs, 0, NULL, NULL);
    if (cl_err)
        err(cl_error_string(cl_err), userData);
    clEnqueueReadBuffer(commands, buff, CL_TRUE, 0, bytes, out, 0, NULL, NULL);
    return cl_err;
    }

    int Cldconv::convolution(float *out, float *in1, float *in2) {
    size_t bytes = vsize * sizeof(cl_float);
    if (wp > irsize) {
        int front = wp - irsize;
        bytes = (vsize - front) * sizeof(cl_float);
        clEnqueueWriteBuffer(commands, coefs, CL_TRUE, wp * sizeof(cl_float), bytes,
                            in2, 0, NULL, NULL);
        bytes = front * sizeof(cl_float);
        clEnqueueWriteBuffer(commands, coefs, CL_TRUE, 0, bytes,
                            &in2[vsize - front], 0, NULL, NULL);
    } else
        clEnqueueWriteBuffer(commands, coefs, CL_TRUE, wp * sizeof(cl_float), bytes,
                            in2, 0, NULL, NULL);
    return convolution(out, in1);
    }

    int Cldconv::push_ir(float *ir) {
    return clEnqueueWriteBuffer(commands, coefs, CL_TRUE, 0,
                                irsize * sizeof(float), ir, 0, NULL, NULL);
    }
}
