double Cheby(int n, double x)
{
    if (n == 0)
    {
        return 1;
    }
    else if (n == 1)
    {
        return x;
    }
    else
    {
        return (2 * x * Cheby(n - 1, x) - Cheby(n - 2, x));
    }
}