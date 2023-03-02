%{

%}

%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(const std::runtime_error& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(...)
  {
    std::cout << "unknown exception\n";
  }
}
