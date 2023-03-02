bool comparator(double x, double y)
{
    if( x > y) return true;
    return false;
}

float opamp_comparator(double x, double y, float Vm=-1.0f, float Vp=1.0f)
{
    if(x < y) return Vm;
    return Vp;
}

// opamp comparator
// if(last_x < 0 && x > 0) crossed
// if(last_x > 0 && x < 0) crossed
bool ZeroCrossingDetection(float x) {
    return x == 0.0f;
}
