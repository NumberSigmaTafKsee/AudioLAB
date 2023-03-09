#pragma once

namespace Casino::IPP
{
    inline int kaiserTapsEstimate(double delta, double As)
    {
        return int((As - 8.0) / (2.285 * delta) + 0.5);
    }
    struct Resample32
    {
        IppsResamplingPolyphaseFixed_32f * pSpec;
        size_t N;
        Ipp64f Time;
        int Outlen;
        Resample32(int inRate, int outRate, Ipp32f rollf=0.9f, Ipp32f as = 80.0f) {
            int size,height,length;
            Ipp32f alpha = kaiserBeta(as);
            Ipp32f delta = (1.0f - rollf) * M_PI;
            int n = kaiserTapsEstimate(delta,as);
            IppStatus status = ippsResamplePolyphaseFixedGetSize_32f(inRate,outRate,n,&size,&length,&height,ippAlgHintFast);
            checkStatus(status);
            pSpec = (IppsResamplingPolyphaseFixed_32f*)ippsMalloc_8u(size);
            status = ippsResamplePolyphaseFixedInit_32f(inRate,outRate,n,rollf,alpha,(IppsResamplingPolyphaseFixed_32f*)pSpec,ippAlgHintFast);
            checkStatus(status);            
        }
        ~Resample32() {
            if(pSpec) ippsFree(pSpec);
        }
        void Execute(const Ipp32f * pSrc, int len, Ipp32f* pDst, Ipp64f norm=0.98) {
            IppStatus status = ippsResamplePolyphaseFixed_32f(pSrc,len,pDst,norm,&Time,&Outlen,(IppsResamplingPolyphaseFixed_32f*)pSpec);
            checkStatus(status);
        }
    };
}