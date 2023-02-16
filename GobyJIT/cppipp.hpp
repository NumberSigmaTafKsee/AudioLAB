#pragma once

#include "ipp.h"
#include "ippe.h"

namespace ipp_cpp {
	inline namespace ippch {
		inline auto RegExpFind(const Ipp8u* pSrc, int srcLen, IppRegExpState* pRegExpState, IppRegExpFind* pFind, int* pNumFind) {
			return ippsRegExpFind_8u(pSrc, srcLen, pRegExpState, pFind, pNumFind);
		}

		inline auto ConvertUTF(const Ipp8u* pSrc, Ipp32u* pSrcLen, Ipp16u* pDst, Ipp32u* pDstLen, int BEFlag) {
			return ippsConvertUTF_8u16u(pSrc, pSrcLen, pDst, pDstLen, BEFlag);
		}

		inline auto ConvertUTF(const Ipp16u* pSrc, Ipp32u* pSrcLen, Ipp8u* pDst, Ipp32u* pDstLen, int BEFlag) {
			return ippsConvertUTF_16u8u(pSrc, pSrcLen, pDst, pDstLen, BEFlag);
		}

		inline auto RegExpReplace(const Ipp8u* pSrc, int* pSrcLenOffset, Ipp8u* pDst, int* pDstLen, IppRegExpFind* pFind, int* pNumFind, IppRegExpState* pRegExpState, IppRegExpReplaceState* pReplaceState) {
			return ippsRegExpReplace_8u(pSrc, pSrcLenOffset, pDst, pDstLen, pFind, pNumFind, pRegExpState, pReplaceState);
		}
	}

	inline namespace ippdc {
		inline auto MTFInit(IppMTFState_8u* pMTFState) {
			return ippsMTFInit_8u(pMTFState);
		}

		template <typename T> 
		auto MTFGetSize(int* pMTFStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto MTFGetSize<Ipp8u>(int* pMTFStateSize) {
			return ippsMTFGetSize_8u(pMTFStateSize);
		}

		inline auto MTFFwd(const Ipp8u* pSrc, Ipp8u* pDst, int len, IppMTFState_8u* pMTFState) {
			return ippsMTFFwd_8u(pSrc, pDst, len, pMTFState);
		}

		inline auto MTFInv(const Ipp8u* pSrc, Ipp8u* pDst, int len, IppMTFState_8u* pMTFState) {
			return ippsMTFInv_8u(pSrc, pDst, len, pMTFState);
		}

		template <typename T> auto BWTFwdGetSize(int wndSize, int* pBWTFwdBuffSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto BWTFwdGetSize<Ipp8u>(int wndSize, int* pBWTFwdBuffSize) {
			return ippsBWTFwdGetSize_8u(wndSize, pBWTFwdBuffSize);
		}

		inline auto BWTFwd(const Ipp8u* pSrc, Ipp8u* pDst, int len, int* index, Ipp8u* pBWTFwdBuff) {
			return ippsBWTFwd_8u(pSrc, pDst, len, index, pBWTFwdBuff);
		}

		template <typename T> 
		auto BWTFwdGetBufSize_SelectSort(Ipp32u wndSize, Ipp32u* pBWTFwdBufSize, IppBWTSortAlgorithmHint sortAlgorithmHint) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto BWTFwdGetBufSize_SelectSort<Ipp8u>(Ipp32u wndSize, Ipp32u* pBWTFwdBufSize, IppBWTSortAlgorithmHint sortAlgorithmHint) {
			return ippsBWTFwdGetBufSize_SelectSort_8u(wndSize, pBWTFwdBufSize, sortAlgorithmHint);
		}

		inline auto BWTFwd_SelectSort(const Ipp8u* pSrc, Ipp8u* pDst, Ipp32u len, Ipp32u* index, Ipp8u* pBWTFwdBuf, IppBWTSortAlgorithmHint sortAlgorithmHint) {
			return ippsBWTFwd_SelectSort_8u(pSrc, pDst, len, index, pBWTFwdBuf, sortAlgorithmHint);
		}

		template <typename T> auto BWTInvGetSize(int wndSize, int* pBWTInvBuffSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto BWTInvGetSize<Ipp8u>(int wndSize, int* pBWTInvBuffSize) {
			return ippsBWTInvGetSize_8u(wndSize, pBWTInvBuffSize);
		}

		inline auto BWTInv(const Ipp8u* pSrc, Ipp8u* pDst, int len, int index, Ipp8u* pBWTInvBuff) {
			return ippsBWTInv_8u(pSrc, pDst, len, index, pBWTInvBuff);
		}

		template <typename T> auto LZSSGetSize(int* pLZSSStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto LZSSGetSize<Ipp8u>(int* pLZSSStateSize) {
			return ippsLZSSGetSize_8u(pLZSSStateSize);
		}

		inline auto EncodeLZSSInit(IppLZSSState_8u* pLZSSState) {
			return ippsEncodeLZSSInit_8u(pLZSSState);
		}

		inline auto EncodeLZSS(Ipp8u** ppSrc, int* pSrcLen, Ipp8u** ppDst, int* pDstLen, IppLZSSState_8u* pLZSSState) {
			return ippsEncodeLZSS_8u(ppSrc, pSrcLen, ppDst, pDstLen, pLZSSState);
		}

		inline auto EncodeLZSSFlush(Ipp8u** ppDst, int* pDstLen, IppLZSSState_8u* pLZSSState) {
			return ippsEncodeLZSSFlush_8u(ppDst, pDstLen, pLZSSState);
		}

		inline auto DecodeLZSSInit(IppLZSSState_8u* pLZSSState) {
			return ippsDecodeLZSSInit_8u(pLZSSState);
		}

		inline auto DecodeLZSS(Ipp8u** ppSrc, int* pSrcLen, Ipp8u** ppDst, int* pDstLen, IppLZSSState_8u* pLZSSState) {
			return ippsDecodeLZSS_8u(ppSrc, pSrcLen, ppDst, pDstLen, pLZSSState);
		}

		inline auto Inflate(Ipp8u** ppSrc, unsigned int* pSrcLen, Ipp32u* pCode, unsigned int* pCodeLenBits, unsigned int winIdx, Ipp8u** ppDst, unsigned int* pDstLen, unsigned int dstIdx, IppInflateMode* pMode, IppInflateState* pIppInflateState) {
			return ippsInflate_8u(ppSrc, pSrcLen, pCode, pCodeLenBits, winIdx, ppDst, pDstLen, dstIdx, pMode, pIppInflateState);
		}

		template <typename T> auto RLEGetSize_BZ2(int* pRLEStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto RLEGetSize_BZ2<Ipp8u>(int* pRLEStateSize) {
			return ippsRLEGetSize_BZ2_8u(pRLEStateSize);
		}

		template <typename T> auto EncodeRLEInit_BZ2(IppRLEState_BZ2* pRLEState) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto EncodeRLEInit_BZ2<Ipp8u>(IppRLEState_BZ2* pRLEState) {
			return ippsEncodeRLEInit_BZ2_8u(pRLEState);
		}

		inline auto EncodeRLE_BZ2(Ipp8u** ppSrc, int* pSrcLen, Ipp8u* pDst, int* pDstLen, IppRLEState_BZ2* pRLEState) {
			return ippsEncodeRLE_BZ2_8u(ppSrc, pSrcLen, pDst, pDstLen, pRLEState);
		}

		inline auto EncodeRLEFlush_BZ2(Ipp8u* pDst, int* pDstLen, IppRLEState_BZ2* pRLEState) {
			return ippsEncodeRLEFlush_BZ2_8u(pDst, pDstLen, pRLEState);
		}

		template <typename T> auto DecodeRLEStateInit_BZ2(IppRLEState_BZ2* pRLEState) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DecodeRLEStateInit_BZ2<Ipp8u>(IppRLEState_BZ2* pRLEState) {
			return ippsDecodeRLEStateInit_BZ2_8u(pRLEState);
		}

		inline auto DecodeRLEState_BZ2(Ipp8u** ppSrc, Ipp32u* pSrcLen, Ipp8u** ppDst, Ipp32u* pDstLen, IppRLEState_BZ2* pRLEState) {
			return ippsDecodeRLEState_BZ2_8u(ppSrc, pSrcLen, ppDst, pDstLen, pRLEState);
		}

		inline auto DecodeRLEStateFlush_BZ2(IppRLEState_BZ2* pRLEState, Ipp8u** ppDst, Ipp32u* pDstLen) {
			return ippsDecodeRLEStateFlush_BZ2_8u(pRLEState, ppDst, pDstLen);
		}

		inline auto RLEGetInUseTable(Ipp8u inUse[256], IppRLEState_BZ2* pRLEState) {
			return ippsRLEGetInUseTable_8u(inUse, pRLEState);
		}

		inline auto EncodeZ1Z2_BZ2(Ipp8u** ppSrc, int* pSrcLen, Ipp16u* pDst, int* pDstLen, int freqTable[258]) {
			return ippsEncodeZ1Z2_BZ2_8u16u(ppSrc, pSrcLen, pDst, pDstLen, freqTable);
		}

		inline auto DecodeZ1Z2_BZ2(Ipp16u** ppSrc, int* pSrcLen, Ipp8u* pDst, int* pDstLen) {
			return ippsDecodeZ1Z2_BZ2_16u8u(ppSrc, pSrcLen, pDst, pDstLen);
		}

		inline auto ReduceDictionary_8u(const Ipp8u inUse[256], Ipp8u* pSrcDst, int srcDstLen, int* pSizeDictionary) {
			return ippsReduceDictionary_8u_I(inUse, pSrcDst, srcDstLen, pSizeDictionary);
		}

		inline auto ExpandDictionary_8u(const Ipp8u inUse[256], Ipp8u* pSrcDst, int srcDstLen, int sizeDictionary) {
			return ippsExpandDictionary_8u_I(inUse, pSrcDst, srcDstLen, sizeDictionary);
		}

		inline auto CRC32_BZ2(const Ipp8u* pSrc, int srcLen, Ipp32u* pCRC32) {
			return ippsCRC32_BZ2_8u(pSrc, srcLen, pCRC32);
		}

		template <typename T> auto EncodeHuffGetSize_BZ2(int wndSize, int* pEncodeHuffStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto EncodeHuffGetSize_BZ2<Ipp16u>(int wndSize, int* pEncodeHuffStateSize) {
			return ippsEncodeHuffGetSize_BZ2_16u8u(wndSize, pEncodeHuffStateSize);
		}

		inline auto EncodeHuffInit_BZ2(int sizeDictionary, const int freqTable[258], const Ipp16u* pSrc, int srcLen, IppEncodeHuffState_BZ2* pEncodeHuffState) {
			return ippsEncodeHuffInit_BZ2_16u8u(sizeDictionary, freqTable, pSrc, srcLen, pEncodeHuffState);
		}

		template <typename T> auto PackHuffContext_BZ2(Ipp32u* pCode, int* pCodeLenBits, Ipp8u* pDst, int* pDstLen, IppEncodeHuffState_BZ2* pEncodeHuffState) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto PackHuffContext_BZ2<Ipp16u>(Ipp32u* pCode, int* pCodeLenBits, Ipp8u* pDst, int* pDstLen, IppEncodeHuffState_BZ2* pEncodeHuffState) {
			return ippsPackHuffContext_BZ2_16u8u(pCode, pCodeLenBits, pDst, pDstLen, pEncodeHuffState);
		}

		inline auto EncodeHuff_BZ2(Ipp32u* pCode, int* pCodeLenBits, Ipp16u** ppSrc, int* pSrcLen, Ipp8u* pDst, int* pDstLen, IppEncodeHuffState_BZ2* pEncodeHuffState) {
			return ippsEncodeHuff_BZ2_16u8u(pCode, pCodeLenBits, ppSrc, pSrcLen, pDst, pDstLen, pEncodeHuffState);
		}

		template <typename T> auto DecodeHuffGetSize_BZ2(int wndSize, int* pDecodeHuffStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DecodeHuffGetSize_BZ2<Ipp8u>(int wndSize, int* pDecodeHuffStateSize) {
			return ippsDecodeHuffGetSize_BZ2_8u16u(wndSize, pDecodeHuffStateSize);
		}

		template <typename T> auto DecodeHuffInit_BZ2(int sizeDictionary, IppDecodeHuffState_BZ2* pDecodeHuffState) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DecodeHuffInit_BZ2<Ipp8u>(int sizeDictionary, IppDecodeHuffState_BZ2* pDecodeHuffState) {
			return ippsDecodeHuffInit_BZ2_8u16u(sizeDictionary, pDecodeHuffState);
		}

		inline auto UnpackHuffContext_BZ2(Ipp32u* pCode, int* pCodeLenBits, Ipp8u** ppSrc, int* pSrcLen, IppDecodeHuffState_BZ2* pDecodeHuffState) {
			return ippsUnpackHuffContext_BZ2_8u16u(pCode, pCodeLenBits, ppSrc, pSrcLen, pDecodeHuffState);
		}

		inline auto DecodeHuff_BZ2(Ipp32u* pCode, int* pCodeLenBits, Ipp8u** ppSrc, int* pSrcLen, Ipp16u* pDst, int* pDstLen, IppDecodeHuffState_BZ2* pDecodeHuffState) {
			return ippsDecodeHuff_BZ2_8u16u(pCode, pCodeLenBits, ppSrc, pSrcLen, pDst, pDstLen, pDecodeHuffState);
		}

		template <typename T> auto DecodeBlockGetSize_BZ2(int blockSize, int* pBuffSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DecodeBlockGetSize_BZ2<Ipp8u>(int blockSize, int* pBuffSize) {
			return ippsDecodeBlockGetSize_BZ2_8u(blockSize, pBuffSize);
		}

		inline auto DecodeBlock_BZ2(const Ipp16u* pSrc, int srcLen, Ipp8u* pDst, int* pDstLen, int index, int dictSize, const Ipp8u inUse[256], Ipp8u* pBuff) {
			return ippsDecodeBlock_BZ2_16u8u(pSrc, srcLen, pDst, pDstLen, index, dictSize, inUse, pBuff);
		}

		inline auto EncodeLZO(const Ipp8u* pSrc, Ipp32u srcLen, Ipp8u* pDst, Ipp32u* pDstLen, IppLZOState_8u* pLZOState) {
			return ippsEncodeLZO_8u(pSrc, srcLen, pDst, pDstLen, pLZOState);
		}

		inline auto EncodeLZOInit(IppLZOMethod method, Ipp32u maxInputLen, IppLZOState_8u* pLZOState) {
			return ippsEncodeLZOInit_8u(method, maxInputLen, pLZOState);
		}

		inline auto DecodeLZO(const Ipp8u* pSrc, Ipp32u srcLen, Ipp8u* pDst, Ipp32u* pDstLen) {
			return ippsDecodeLZO_8u(pSrc, srcLen, pDst, pDstLen);
		}

		inline auto DecodeLZOSafe(const Ipp8u* pSrc, Ipp32u srcLen, Ipp8u* pDst, Ipp32u* pDstLen) {
			return ippsDecodeLZOSafe_8u(pSrc, srcLen, pDst, pDstLen);
		}

		template <typename T> auto EncodeLZ4HashTableGetSize(int* pHashTableSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto EncodeLZ4HashTableGetSize<Ipp8u>(int* pHashTableSize) {
			return ippsEncodeLZ4HashTableGetSize_8u(pHashTableSize);
		}

		inline auto EncodeLZ4HashTableInit(Ipp8u* pHashTable, int srcLen) {
			return ippsEncodeLZ4HashTableInit_8u(pHashTable, srcLen);
		}

		inline auto EncodeLZ4DictHashTableInit(Ipp8u* pHashTable, int srcLen) {
			return ippsEncodeLZ4DictHashTableInit_8u(pHashTable, srcLen);
		}

		inline auto EncodeLZ4LoadDict(Ipp8u* pHashTable, const Ipp8u* pDict, int dictLen) {
			return ippsEncodeLZ4LoadDict_8u(pHashTable, pDict, dictLen);
		}

		inline auto EncodeLZ4(const Ipp8u* pSrc, int srcLen, Ipp8u* pDst, int* pDstLen, Ipp8u* pHashTable) {
			return ippsEncodeLZ4_8u(pSrc, srcLen, pDst, pDstLen, pHashTable);
		}

		inline auto EncodeLZ4Fast(const Ipp8u* pSrc, int srcLen, Ipp8u* pDst, int* pDstLen, Ipp8u* pHashTable, int acceleration) {
			return ippsEncodeLZ4Fast_8u(pSrc, srcLen, pDst, pDstLen, pHashTable, acceleration);
		}

		inline auto EncodeLZ4Safe(const Ipp8u* pSrc, int* srcLen, Ipp8u* pDst, int* pDstLen, Ipp8u* pHashTable) {
			return ippsEncodeLZ4Safe_8u(pSrc, srcLen, pDst, pDstLen, pHashTable);
		}

		inline auto EncodeLZ4Dict(const Ipp8u* pSrc, int srcIdx, int srcLen, Ipp8u* pDst, int* pDstLen, Ipp8u* pHashTable, const Ipp8u* pDict, int dictLen) {
			return ippsEncodeLZ4Dict_8u(pSrc, srcIdx, srcLen, pDst, pDstLen, pHashTable, pDict, dictLen);
		}

		inline auto EncodeLZ4DictSafe(const Ipp8u* pSrc, int srcIdx, int* pSrcLen, Ipp8u* pDst, int* pDstLen, Ipp8u* pHashTable, const Ipp8u* pDict, int dictLen) {
			return ippsEncodeLZ4DictSafe_8u(pSrc, srcIdx, pSrcLen, pDst, pDstLen, pHashTable, pDict, dictLen);
		}

		template <typename T> auto EncodeLZ4HCHashTableGetSize(int* pHashTableSize, int* pPrevTableSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto EncodeLZ4HCHashTableGetSize<Ipp8u>(int* pHashTableSize, int* pPrevTableSize) {
			return ippsEncodeLZ4HCHashTableGetSize_8u(pHashTableSize, pPrevTableSize);
		}

		inline auto EncodeLZ4HCHashTableInit(Ipp8u** ppHashTables) {
			return ippsEncodeLZ4HCHashTableInit_8u(ppHashTables);
		}

		inline auto EncodeLZ4HC(const Ipp8u* pSrc, int srcIdx, int* pSrcLen, Ipp8u* pDst, int* pDstLen, Ipp8u** ppHashTables, const Ipp8u* pDict, int dictLen, int level) {
			return ippsEncodeLZ4HC_8u(pSrc, srcIdx, pSrcLen, pDst, pDstLen, ppHashTables, pDict, dictLen, level);
		}

		inline auto EncodeLZ4HCDictLimit(const Ipp8u* pSrc, int srcIdx, int* pSrcLen, Ipp8u* pDst, int* pDstLen, Ipp8u** ppHashTables, const Ipp8u* pDict, int dictLen, int level, int lowDictIdx) {
			return ippsEncodeLZ4HCDictLimit_8u(pSrc, srcIdx, pSrcLen, pDst, pDstLen, ppHashTables, pDict, dictLen, level, lowDictIdx);
		}

		inline auto DecodeLZ4(const Ipp8u* pSrc, int srcLen, Ipp8u* pDst, int* pDstLen) {
			return ippsDecodeLZ4_8u(pSrc, srcLen, pDst, pDstLen);
		}

		inline auto DecodeLZ4Dict(const Ipp8u* pSrc, int* pSrcLen, Ipp8u* pDst, int dstIdx, int* pDstLen, const Ipp8u* pDict, int dictSize) {
			return ippsDecodeLZ4Dict_8u(pSrc, pSrcLen, pDst, dstIdx, pDstLen, pDict, dictSize);
		}

		template <typename T> auto EncodeZfpGetStateSize(int* pStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto EncodeZfpGetStateSize<Ipp32f>(int* pStateSize) {
			return ippsEncodeZfpGetStateSize_32f(pStateSize);
		}

		inline auto EncodeZfpInit(Ipp8u* pDst, int dstLen, IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfpInit_32f(pDst, dstLen, pState);
		}

		inline auto EncodeZfpInitLong(Ipp8u* pDst, Ipp64u dstLen, IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfpInitLong_32f(pDst, dstLen, pState);
		}

		inline auto EncodeZfpSet(int minBits, int maxBits, int maxPrec, int minExp, IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfpSet_32f(minBits, maxBits, maxPrec, minExp, pState);
		}

		inline auto EncodeZfpSetAccuracy(Ipp64f precision, IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfpSetAccuracy_32f(precision, pState);
		}

		inline auto EncodeZfp444(const Ipp32f* pSrc, int srcStep, int srcPlaneStep, IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfp444_32f(pSrc, srcStep, srcPlaneStep, pState);
		}

		inline auto EncodeZfpFlush(IppEncodeZfpState_32f* pState) {
			return ippsEncodeZfpFlush_32f(pState);
		}

		inline auto EncodeZfpGetCompressedSize(IppEncodeZfpState_32f* pState, int* pCompressedSize) {
			return ippsEncodeZfpGetCompressedSize_32f(pState, pCompressedSize);
		}

		inline auto EncodeZfpGetCompressedSizeLong(IppEncodeZfpState_32f* pState, Ipp64u* pCompressedSize) {
			return ippsEncodeZfpGetCompressedSizeLong_32f(pState, pCompressedSize);
		}

		inline auto EncodeZfpGetCompressedBitSize(IppEncodeZfpState_32f* pState, Ipp64u* pCompressedBitSize) {
			return ippsEncodeZfpGetCompressedBitSize_32f(pState, pCompressedBitSize);
		}

		template <typename T> auto DecodeZfpGetStateSize(int* pStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DecodeZfpGetStateSize<Ipp32f>(int* pStateSize) {
			return ippsDecodeZfpGetStateSize_32f(pStateSize);
		}

		inline auto DecodeZfpInit(const Ipp8u* pSrc, int srcLen, IppDecodeZfpState_32f* pState) {
			return ippsDecodeZfpInit_32f(pSrc, srcLen, pState);
		}

		inline auto DecodeZfpInitLong(const Ipp8u* pSrc, Ipp64u srcLen, IppDecodeZfpState_32f* pState) {
			return ippsDecodeZfpInitLong_32f(pSrc, srcLen, pState);
		}

		inline auto DecodeZfpSetAccuracy(Ipp64f precision, IppDecodeZfpState_32f* pState) {
			return ippsDecodeZfpSetAccuracy_32f(precision, pState);
		}

		inline auto DecodeZfpSet(int minBits, int maxBits, int maxPrec, int minExp, IppDecodeZfpState_32f* pState) {
			return ippsDecodeZfpSet_32f(minBits, maxBits, maxPrec, minExp, pState);
		}

		inline auto DecodeZfp444(IppDecodeZfpState_32f* pState, Ipp32f* pDst, int dstStep, int dstPlaneStep) {
			return ippsDecodeZfp444_32f(pState, pDst, dstStep, dstPlaneStep);
		}

		inline auto DecodeZfpGetDecompressedSize(IppDecodeZfpState_32f* pState, int* pDecompressedSize) {
			return ippsDecodeZfpGetDecompressedSize_32f(pState, pDecompressedSize);
		}

		inline auto DecodeZfpGetDecompressedSizeLong(IppDecodeZfpState_32f* pState, Ipp64u* pDecompressedSize) {
			return ippsDecodeZfpGetDecompressedSizeLong_32f(pState, pDecompressedSize);
		}
	}

	inline namespace ippe {
		inline auto MimoMMSE_1X2(Ipp16sc* pSrcH[2], int srcHStride2, int srcHStride1, int srcHStride0, Ipp16sc* pSrcY[4][12], int Sigma2, IppFourSymb* pDstX, int dstXStride1, int dstXStride0, int numSymb, int numSC, int SINRIdx, Ipp32f* pDstSINR, int scaleFactor) {
			return ippsMimoMMSE_1X2_16sc(pSrcH, srcHStride2, srcHStride1, srcHStride0, pSrcY, Sigma2, pDstX, dstXStride1, dstXStride0, numSymb, numSC, SINRIdx, pDstSINR, scaleFactor);
		}

		inline auto MimoMMSE_2X2(Ipp16sc* pSrcH[2], int srcHStride2, int srcHStride1, int srcHStride0, Ipp16sc* pSrcY[4][12], int Sigma2, IppFourSymb* pDstX, int dstXStride1, int dstXStride0, int numSymb, int numSC, int SINRIdx, Ipp32f* pDstSINR, int scaleFactor) {
			return ippsMimoMMSE_2X2_16sc(pSrcH, srcHStride2, srcHStride1, srcHStride0, pSrcY, Sigma2, pDstX, dstXStride1, dstXStride0, numSymb, numSC, SINRIdx, pDstSINR, scaleFactor);
		}

		inline auto MimoMMSE_1X4(Ipp16sc* pSrcH[2], int srcHStride2, int srcHStride1, int srcHStride0, Ipp16sc* pSrcY[4][12], int Sigma2, IppFourSymb* pDstX, int dstXStride1, int dstXStride0, int numSymb, int numSC, int SINRIdx, Ipp32f* pDstSINR, int scaleFactor) {
			return ippsMimoMMSE_1X4_16sc(pSrcH, srcHStride2, srcHStride1, srcHStride0, pSrcY, Sigma2, pDstX, dstXStride1, dstXStride0, numSymb, numSC, SINRIdx, pDstSINR, scaleFactor);
		}

		inline auto MimoMMSE_2X4(Ipp16sc* pSrcH[2], int srcHStride2, int srcHStride1, int srcHStride0, Ipp16sc* pSrcY[4][12], int Sigma2, IppFourSymb* pDstX, int dstXStride1, int dstXStride0, int numSymb, int numSC, int SINRIdx, Ipp32f* pDstSINR, int scaleFactor) {
			return ippsMimoMMSE_2X4_16sc(pSrcH, srcHStride2, srcHStride1, srcHStride0, pSrcY, Sigma2, pDstX, dstXStride1, dstXStride0, numSymb, numSC, SINRIdx, pDstSINR, scaleFactor);
		}

		inline auto CRC24a(Ipp8u* pSrc, int len, Ipp32u* pCRC24) {
			return ippsCRC24a_8u(pSrc, len, pCRC24);
		}

		inline auto CRC24b(Ipp8u* pSrc, int len, Ipp32u* pCRC24) {
			return ippsCRC24b_8u(pSrc, len, pCRC24);
		}

		inline auto CRC24c(Ipp8u* pSrc, int len, Ipp32u* pCRC24) {
			return ippsCRC24c_8u(pSrc, len, pCRC24);
		}

		inline auto CRC16(Ipp8u* pSrc, int len, Ipp32u* pCRC16) {
			return ippsCRC16_8u(pSrc, len, pCRC16);
		}

		inline auto CRC24a(Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int bitLen) {
			return ippsCRC24a_1u(pSrc, srcBitOffset, pDst, dstBitOffset, bitLen);
		}

		inline auto CRC24b(Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int bitLen) {
			return ippsCRC24b_1u(pSrc, srcBitOffset, pDst, dstBitOffset, bitLen);
		}

		inline auto CRC24c(Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int bitLen) {
			return ippsCRC24c_1u(pSrc, srcBitOffset, pDst, dstBitOffset, bitLen);
		}

		inline auto CRC16(Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int bitLen) {
			return ippsCRC16_1u(pSrc, srcBitOffset, pDst, dstBitOffset, bitLen);
		}

		inline auto CRC(const Ipp8u* pSrc, int len, Ipp64u poly, const Ipp8u optPoly[128], Ipp32u init, Ipp32u* pCRC) {
			return ippsCRC_8u(pSrc, len, poly, optPoly, init, pCRC);
		}

		inline auto GenCRCOptPoly(Ipp64u poly, Ipp8u optPoly[128]) {
			return ippsGenCRCOptPoly_8u(poly, optPoly);
		}
	}

	inline namespace ipps {
		struct Ipp24u {
			Ipp8u data[3];
		};

		struct Ipp24s {
			Ipp8u data[3];
		};

		// todo: implement proper fp16
		struct Ipp16f {
			Ipp16s data;
		};

		inline auto Copy(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
			return ippsCopy_8u(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsCopy_16s(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsCopy_16sc(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsCopy_32f(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsCopy_32fc(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsCopy_64f(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsCopy_64fc(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
			return ippsCopy_32s(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
			return ippsCopy_32sc(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
			return ippsCopy_64s(pSrc, pDst, len);
		}

		inline auto Copy(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
			return ippsCopy_64sc(pSrc, pDst, len);
		}

		inline auto CopyLE(const Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int len) {
			return ippsCopyLE_1u(pSrc, srcBitOffset, pDst, dstBitOffset, len);
		}

		inline auto CopyBE(const Ipp8u* pSrc, int srcBitOffset, Ipp8u* pDst, int dstBitOffset, int len) {
			return ippsCopyBE_1u(pSrc, srcBitOffset, pDst, dstBitOffset, len);
		}

		inline auto Move(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
			return ippsMove_8u(pSrc, pDst, len);
		}

		inline auto Move(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsMove_16s(pSrc, pDst, len);
		}

		inline auto Move(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsMove_16sc(pSrc, pDst, len);
		}

		inline auto Move(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsMove_32f(pSrc, pDst, len);
		}

		inline auto Move(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsMove_32fc(pSrc, pDst, len);
		}

		inline auto Move(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsMove_64f(pSrc, pDst, len);
		}

		inline auto Move(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsMove_64fc(pSrc, pDst, len);
		}

		inline auto Move(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
			return ippsMove_32s(pSrc, pDst, len);
		}

		inline auto Move(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
			return ippsMove_32sc(pSrc, pDst, len);
		}

		inline auto Move(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
			return ippsMove_64s(pSrc, pDst, len);
		}

		inline auto Move(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
			return ippsMove_64sc(pSrc, pDst, len);
		}

		inline auto Set(Ipp8u val, Ipp8u* pDst, int len) {
			return ippsSet_8u(val, pDst, len);
		}

		inline auto Set(Ipp16s val, Ipp16s* pDst, int len) {
			return ippsSet_16s(val, pDst, len);
		}

		inline auto Set(Ipp16sc val, Ipp16sc* pDst, int len) {
			return ippsSet_16sc(val, pDst, len);
		}

		inline auto Set(Ipp32s val, Ipp32s* pDst, int len) {
			return ippsSet_32s(val, pDst, len);
		}

		inline auto Set(Ipp32sc val, Ipp32sc* pDst, int len) {
			return ippsSet_32sc(val, pDst, len);
		}

		inline auto Set(Ipp32f val, Ipp32f* pDst, int len) {
			return ippsSet_32f(val, pDst, len);
		}

		inline auto Set(Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsSet_32fc(val, pDst, len);
		}

		inline auto Set(Ipp64s val, Ipp64s* pDst, int len) {
			return ippsSet_64s(val, pDst, len);
		}

		inline auto Set(Ipp64sc val, Ipp64sc* pDst, int len) {
			return ippsSet_64sc(val, pDst, len);
		}

		inline auto Set(Ipp64f val, Ipp64f* pDst, int len) {
			return ippsSet_64f(val, pDst, len);
		}

		inline auto Set(Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsSet_64fc(val, pDst, len);
		}

		inline auto Zero(Ipp8u* pDst, int len) {
			return ippsZero_8u(pDst, len);
		}

		inline auto Zero(Ipp16s* pDst, int len) {
			return ippsZero_16s(pDst, len);
		}

		inline auto Zero(Ipp16sc* pDst, int len) {
			return ippsZero_16sc(pDst, len);
		}

		inline auto Zero(Ipp32f* pDst, int len) {
			return ippsZero_32f(pDst, len);
		}

		inline auto Zero(Ipp32fc* pDst, int len) {
			return ippsZero_32fc(pDst, len);
		}

		inline auto Zero(Ipp64f* pDst, int len) {
			return ippsZero_64f(pDst, len);
		}

		inline auto Zero(Ipp64fc* pDst, int len) {
			return ippsZero_64fc(pDst, len);
		}

		inline auto Zero(Ipp32s* pDst, int len) {
			return ippsZero_32s(pDst, len);
		}

		inline auto Zero(Ipp32sc* pDst, int len) {
			return ippsZero_32sc(pDst, len);
		}

		inline auto Zero(Ipp64s* pDst, int len) {
			return ippsZero_64s(pDst, len);
		}

		inline auto Zero(Ipp64sc* pDst, int len) {
			return ippsZero_64sc(pDst, len);
		}

		inline auto Tone(Ipp16s* pDst, int len, Ipp16s magn, Ipp32f rFreq, Ipp32f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_16s(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Tone(Ipp16sc* pDst, int len, Ipp16s magn, Ipp32f rFreq, Ipp32f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_16sc(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Tone(Ipp32f* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_32f(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Tone(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_32fc(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Tone(Ipp64f* pDst, int len, Ipp64f magn, Ipp64f rFreq, Ipp64f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_64f(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Tone(Ipp64fc* pDst, int len, Ipp64f magn, Ipp64f rFreq, Ipp64f* pPhase, IppHintAlgorithm hint) {
			return ippsTone_64fc(pDst, len, magn, rFreq, pPhase, hint);
		}

		inline auto Triangle(Ipp64f* pDst, int len, Ipp64f magn, Ipp64f rFreq, Ipp64f asym, Ipp64f* pPhase) {
			return ippsTriangle_64f(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto Triangle(Ipp64fc* pDst, int len, Ipp64f magn, Ipp64f rFreq, Ipp64f asym, Ipp64f* pPhase) {
			return ippsTriangle_64fc(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto Triangle(Ipp32f* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f asym, Ipp32f* pPhase) {
			return ippsTriangle_32f(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto Triangle(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f asym, Ipp32f* pPhase) {
			return ippsTriangle_32fc(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto Triangle(Ipp16s* pDst, int len, Ipp16s magn, Ipp32f rFreq, Ipp32f asym, Ipp32f* pPhase) {
			return ippsTriangle_16s(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto Triangle(Ipp16sc* pDst, int len, Ipp16s magn, Ipp32f rFreq, Ipp32f asym, Ipp32f* pPhase) {
			return ippsTriangle_16sc(pDst, len, magn, rFreq, asym, pPhase);
		}

		inline auto RandUniform(Ipp8u* pDst, int len, IppsRandUniState_8u* pRandUniState) {
			return ippsRandUniform_8u(pDst, len, pRandUniState);
		}

		inline auto RandUniform(Ipp16s* pDst, int len, IppsRandUniState_16s* pRandUniState) {
			return ippsRandUniform_16s(pDst, len, pRandUniState);
		}

		inline auto RandUniform(Ipp32f* pDst, int len, IppsRandUniState_32f* pRandUniState) {
			return ippsRandUniform_32f(pDst, len, pRandUniState);
		}

		inline auto RandUniform(Ipp64f* pDst, int len, IppsRandUniState_64f* pRandUniState) {
			return ippsRandUniform_64f(pDst, len, pRandUniState);
		}

		inline auto RandGauss(Ipp8u* pDst, int len, IppsRandGaussState_8u* pRandGaussState) {
			return ippsRandGauss_8u(pDst, len, pRandGaussState);
		}

		inline auto RandGauss(Ipp16s* pDst, int len, IppsRandGaussState_16s* pRandGaussState) {
			return ippsRandGauss_16s(pDst, len, pRandGaussState);
		}

		inline auto RandGauss(Ipp32f* pDst, int len, IppsRandGaussState_32f* pRandGaussState) {
			return ippsRandGauss_32f(pDst, len, pRandGaussState);
		}

		inline auto RandGauss(Ipp64f* pDst, int len, IppsRandGaussState_64f* pRandGaussState) {
			return ippsRandGauss_64f(pDst, len, pRandGaussState);
		}

		template <typename T> auto RandGaussGetSize(int* pRandGaussStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto RandGaussGetSize<Ipp8u>(int* pRandGaussStateSize) {
			return ippsRandGaussGetSize_8u(pRandGaussStateSize);
		}

		template <>
		inline auto RandGaussGetSize<Ipp16s>(int* pRandGaussStateSize) {
			return ippsRandGaussGetSize_16s(pRandGaussStateSize);
		}

		template <>
		inline auto RandGaussGetSize<Ipp32f>(int* pRandGaussStateSize) {
			return ippsRandGaussGetSize_32f(pRandGaussStateSize);
		}

		template <>
		inline auto RandGaussGetSize<Ipp64f>(int* pRandGaussStateSize) {
			return ippsRandGaussGetSize_64f(pRandGaussStateSize);
		}

		inline auto RandGaussInit(IppsRandGaussState_8u* pRandGaussState, Ipp8u mean, Ipp8u stdDev, unsigned int seed) {
			return ippsRandGaussInit_8u(pRandGaussState, mean, stdDev, seed);
		}

		inline auto RandGaussInit(IppsRandGaussState_16s* pRandGaussState, Ipp16s mean, Ipp16s stdDev, unsigned int seed) {
			return ippsRandGaussInit_16s(pRandGaussState, mean, stdDev, seed);
		}

		inline auto RandGaussInit(IppsRandGaussState_32f* pRandGaussState, Ipp32f mean, Ipp32f stdDev, unsigned int seed) {
			return ippsRandGaussInit_32f(pRandGaussState, mean, stdDev, seed);
		}

		inline auto RandGaussInit(IppsRandGaussState_64f* pRandGaussState, Ipp64f mean, Ipp64f stdDev, unsigned int seed) {
			return ippsRandGaussInit_64f(pRandGaussState, mean, stdDev, seed);
		}

		template <typename T> auto RandUniformGetSize(int* pRandUniformStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto RandUniformGetSize<Ipp8u>(int* pRandUniformStateSize) {
			return ippsRandUniformGetSize_8u(pRandUniformStateSize);
		}

		template <>
		inline auto RandUniformGetSize<Ipp16s>(int* pRandUniformStateSize) {
			return ippsRandUniformGetSize_16s(pRandUniformStateSize);
		}

		template <>
		inline auto RandUniformGetSize<Ipp32f>(int* pRandUniformStateSize) {
			return ippsRandUniformGetSize_32f(pRandUniformStateSize);
		}

		template <>
		inline auto RandUniformGetSize<Ipp64f>(int* pRandUniformStateSize) {
			return ippsRandUniformGetSize_64f(pRandUniformStateSize);
		}

		inline auto RandUniformInit(IppsRandUniState_8u* pRandUniState, Ipp8u low, Ipp8u high, unsigned int seed) {
			return ippsRandUniformInit_8u(pRandUniState, low, high, seed);
		}

		inline auto RandUniformInit(IppsRandUniState_16s* pRandUniState, Ipp16s low, Ipp16s high, unsigned int seed) {
			return ippsRandUniformInit_16s(pRandUniState, low, high, seed);
		}

		inline auto RandUniformInit(IppsRandUniState_32f* pRandUniState, Ipp32f low, Ipp32f high, unsigned int seed) {
			return ippsRandUniformInit_32f(pRandUniState, low, high, seed);
		}

		inline auto RandUniformInit(IppsRandUniState_64f* pRandUniState, Ipp64f low, Ipp64f high, unsigned int seed) {
			return ippsRandUniformInit_64f(pRandUniState, low, high, seed);
		}

		inline auto VectorJaehne(Ipp8u* pDst, int len, Ipp8u magn) {
			return ippsVectorJaehne_8u(pDst, len, magn);
		}

		inline auto VectorJaehne(Ipp16u* pDst, int len, Ipp16u magn) {
			return ippsVectorJaehne_16u(pDst, len, magn);
		}

		inline auto VectorJaehne(Ipp16s* pDst, int len, Ipp16s magn) {
			return ippsVectorJaehne_16s(pDst, len, magn);
		}

		inline auto VectorJaehne(Ipp32s* pDst, int len, Ipp32s magn) {
			return ippsVectorJaehne_32s(pDst, len, magn);
		}

		inline auto VectorJaehne(Ipp32f* pDst, int len, Ipp32f magn) {
			return ippsVectorJaehne_32f(pDst, len, magn);
		}

		inline auto VectorJaehne(Ipp64f* pDst, int len, Ipp64f magn) {
			return ippsVectorJaehne_64f(pDst, len, magn);
		}

		inline auto VectorSlope(Ipp8u* pDst, int len, Ipp32f offset, Ipp32f slope) {
			return ippsVectorSlope_8u(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp16u* pDst, int len, Ipp32f offset, Ipp32f slope) {
			return ippsVectorSlope_16u(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp16s* pDst, int len, Ipp32f offset, Ipp32f slope) {
			return ippsVectorSlope_16s(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp32u* pDst, int len, Ipp64f offset, Ipp64f slope) {
			return ippsVectorSlope_32u(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp32s* pDst, int len, Ipp64f offset, Ipp64f slope) {
			return ippsVectorSlope_32s(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp32f* pDst, int len, Ipp32f offset, Ipp32f slope) {
			return ippsVectorSlope_32f(pDst, len, offset, slope);
		}

		inline auto VectorSlope(Ipp64f* pDst, int len, Ipp64f offset, Ipp64f slope) {
			return ippsVectorSlope_64f(pDst, len, offset, slope);
		}

		inline auto AndC_8u(Ipp8u val, Ipp8u* pSrcDst, int len) {
			return ippsAndC_8u_I(val, pSrcDst, len);
		}

		inline auto AndC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len) {
			return ippsAndC_8u(pSrc, val, pDst, len);
		}

		inline auto AndC_16u(Ipp16u val, Ipp16u* pSrcDst, int len) {
			return ippsAndC_16u_I(val, pSrcDst, len);
		}

		inline auto AndC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len) {
			return ippsAndC_16u(pSrc, val, pDst, len);
		}

		inline auto AndC_32u(Ipp32u val, Ipp32u* pSrcDst, int len) {
			return ippsAndC_32u_I(val, pSrcDst, len);
		}

		inline auto AndC(const Ipp32u* pSrc, Ipp32u val, Ipp32u* pDst, int len) {
			return ippsAndC_32u(pSrc, val, pDst, len);
		}

		inline auto And_8u(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len) {
			return ippsAnd_8u_I(pSrc, pSrcDst, len);
		}

		inline auto And(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len) {
			return ippsAnd_8u(pSrc1, pSrc2, pDst, len);
		}

		inline auto And_16u(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len) {
			return ippsAnd_16u_I(pSrc, pSrcDst, len);
		}

		inline auto And(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len) {
			return ippsAnd_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto And_32u(const Ipp32u* pSrc, Ipp32u* pSrcDst, int len) {
			return ippsAnd_32u_I(pSrc, pSrcDst, len);
		}

		inline auto And(const Ipp32u* pSrc1, const Ipp32u* pSrc2, Ipp32u* pDst, int len) {
			return ippsAnd_32u(pSrc1, pSrc2, pDst, len);
		}

		inline auto OrC_8u(Ipp8u val, Ipp8u* pSrcDst, int len) {
			return ippsOrC_8u_I(val, pSrcDst, len);
		}

		inline auto OrC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len) {
			return ippsOrC_8u(pSrc, val, pDst, len);
		}

		inline auto OrC_16u(Ipp16u val, Ipp16u* pSrcDst, int len) {
			return ippsOrC_16u_I(val, pSrcDst, len);
		}

		inline auto OrC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len) {
			return ippsOrC_16u(pSrc, val, pDst, len);
		}

		inline auto OrC_32u(Ipp32u val, Ipp32u* pSrcDst, int len) {
			return ippsOrC_32u_I(val, pSrcDst, len);
		}

		inline auto OrC(const Ipp32u* pSrc, Ipp32u val, Ipp32u* pDst, int len) {
			return ippsOrC_32u(pSrc, val, pDst, len);
		}

		inline auto Or_8u(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len) {
			return ippsOr_8u_I(pSrc, pSrcDst, len);
		}

		inline auto Or(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len) {
			return ippsOr_8u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Or_16u(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len) {
			return ippsOr_16u_I(pSrc, pSrcDst, len);
		}

		inline auto Or(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len) {
			return ippsOr_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Or_32u(const Ipp32u* pSrc, Ipp32u* pSrcDst, int len) {
			return ippsOr_32u_I(pSrc, pSrcDst, len);
		}

		inline auto Or(const Ipp32u* pSrc1, const Ipp32u* pSrc2, Ipp32u* pDst, int len) {
			return ippsOr_32u(pSrc1, pSrc2, pDst, len);
		}

		inline auto XorC_8u(Ipp8u val, Ipp8u* pSrcDst, int len) {
			return ippsXorC_8u_I(val, pSrcDst, len);
		}

		inline auto XorC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len) {
			return ippsXorC_8u(pSrc, val, pDst, len);
		}

		inline auto XorC_16u(Ipp16u val, Ipp16u* pSrcDst, int len) {
			return ippsXorC_16u_I(val, pSrcDst, len);
		}

		inline auto XorC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len) {
			return ippsXorC_16u(pSrc, val, pDst, len);
		}

		inline auto XorC_32u(Ipp32u val, Ipp32u* pSrcDst, int len) {
			return ippsXorC_32u_I(val, pSrcDst, len);
		}

		inline auto XorC(const Ipp32u* pSrc, Ipp32u val, Ipp32u* pDst, int len) {
			return ippsXorC_32u(pSrc, val, pDst, len);
		}

		inline auto Xor_8u(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len) {
			return ippsXor_8u_I(pSrc, pSrcDst, len);
		}

		inline auto Xor(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len) {
			return ippsXor_8u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Xor_16u(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len) {
			return ippsXor_16u_I(pSrc, pSrcDst, len);
		}

		inline auto Xor(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len) {
			return ippsXor_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Xor_32u(const Ipp32u* pSrc, Ipp32u* pSrcDst, int len) {
			return ippsXor_32u_I(pSrc, pSrcDst, len);
		}

		inline auto Xor(const Ipp32u* pSrc1, const Ipp32u* pSrc2, Ipp32u* pDst, int len) {
			return ippsXor_32u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Not_8u(Ipp8u* pSrcDst, int len) {
			return ippsNot_8u_I(pSrcDst, len);
		}

		inline auto Not(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
			return ippsNot_8u(pSrc, pDst, len);
		}

		inline auto Not_16u(Ipp16u* pSrcDst, int len) {
			return ippsNot_16u_I(pSrcDst, len);
		}

		inline auto Not(const Ipp16u* pSrc, Ipp16u* pDst, int len) {
			return ippsNot_16u(pSrc, pDst, len);
		}

		inline auto Not_32u(Ipp32u* pSrcDst, int len) {
			return ippsNot_32u_I(pSrcDst, len);
		}

		inline auto Not(const Ipp32u* pSrc, Ipp32u* pDst, int len) {
			return ippsNot_32u(pSrc, pDst, len);
		}

		inline auto LShiftC_8u(int val, Ipp8u* pSrcDst, int len) {
			return ippsLShiftC_8u_I(val, pSrcDst, len);
		}

		inline auto LShiftC(const Ipp8u* pSrc, int val, Ipp8u* pDst, int len) {
			return ippsLShiftC_8u(pSrc, val, pDst, len);
		}

		inline auto LShiftC_16u(int val, Ipp16u* pSrcDst, int len) {
			return ippsLShiftC_16u_I(val, pSrcDst, len);
		}

		inline auto LShiftC(const Ipp16u* pSrc, int val, Ipp16u* pDst, int len) {
			return ippsLShiftC_16u(pSrc, val, pDst, len);
		}

		inline auto LShiftC_16s(int val, Ipp16s* pSrcDst, int len) {
			return ippsLShiftC_16s_I(val, pSrcDst, len);
		}

		inline auto LShiftC(const Ipp16s* pSrc, int val, Ipp16s* pDst, int len) {
			return ippsLShiftC_16s(pSrc, val, pDst, len);
		}

		inline auto LShiftC_32s(int val, Ipp32s* pSrcDst, int len) {
			return ippsLShiftC_32s_I(val, pSrcDst, len);
		}

		inline auto LShiftC(const Ipp32s* pSrc, int val, Ipp32s* pDst, int len) {
			return ippsLShiftC_32s(pSrc, val, pDst, len);
		}

		inline auto RShiftC_8u(int val, Ipp8u* pSrcDst, int len) {
			return ippsRShiftC_8u_I(val, pSrcDst, len);
		}

		inline auto RShiftC(const Ipp8u* pSrc, int val, Ipp8u* pDst, int len) {
			return ippsRShiftC_8u(pSrc, val, pDst, len);
		}

		inline auto RShiftC_16u(int val, Ipp16u* pSrcDst, int len) {
			return ippsRShiftC_16u_I(val, pSrcDst, len);
		}

		inline auto RShiftC(const Ipp16u* pSrc, int val, Ipp16u* pDst, int len) {
			return ippsRShiftC_16u(pSrc, val, pDst, len);
		}

		inline auto RShiftC_16s(int val, Ipp16s* pSrcDst, int len) {
			return ippsRShiftC_16s_I(val, pSrcDst, len);
		}

		inline auto RShiftC(const Ipp16s* pSrc, int val, Ipp16s* pDst, int len) {
			return ippsRShiftC_16s(pSrc, val, pDst, len);
		}

		inline auto RShiftC_32s(int val, Ipp32s* pSrcDst, int len) {
			return ippsRShiftC_32s_I(val, pSrcDst, len);
		}

		inline auto RShiftC(const Ipp32s* pSrc, int val, Ipp32s* pDst, int len) {
			return ippsRShiftC_32s(pSrc, val, pDst, len);
		}

		inline auto AddC(Ipp8u val, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_8u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsAddC_8u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC_16s(Ipp16s val, Ipp16s* pSrcDst, int len) {
			return ippsAddC_16s_I(val, pSrcDst, len);
		}

		inline auto AddC(Ipp16s val, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_16s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp16s* pSrc, Ipp16s val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsAddC_16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC(Ipp16sc val, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_16sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp16sc* pSrc, Ipp16sc val, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsAddC_16sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC(Ipp16u val, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_16u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsAddC_16u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC(Ipp32s val, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_32s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp32s* pSrc, Ipp32s val, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsAddC_32s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC(Ipp32sc val, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsAddC_32sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp32sc* pSrc, Ipp32sc val, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsAddC_32sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto AddC(const Ipp64u* pSrc, Ipp64u val, Ipp64u* pDst, Ipp32u len, int scaleFactor, IppRoundMode rndMode) {
			return ippsAddC_64u_Sfs(pSrc, val, pDst, len, scaleFactor, rndMode);
		}

		inline auto AddC(const Ipp64s* pSrc, Ipp64s val, Ipp64s* pDst, Ipp32u len, int scaleFactor, IppRoundMode rndMode) {
			return ippsAddC_64s_Sfs(pSrc, val, pDst, len, scaleFactor, rndMode);
		}

		inline auto AddC_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsAddC_32f_I(val, pSrcDst, len);
		}

		inline auto AddC(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsAddC_32f(pSrc, val, pDst, len);
		}

		inline auto AddC_32fc(Ipp32fc val, Ipp32fc* pSrcDst, int len) {
			return ippsAddC_32fc_I(val, pSrcDst, len);
		}

		inline auto AddC(const Ipp32fc* pSrc, Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsAddC_32fc(pSrc, val, pDst, len);
		}

		inline auto AddC_64f(Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsAddC_64f_I(val, pSrcDst, len);
		}

		inline auto AddC(const Ipp64f* pSrc, Ipp64f val, Ipp64f* pDst, int len) {
			return ippsAddC_64f(pSrc, val, pDst, len);
		}

		inline auto AddC_64fc(Ipp64fc val, Ipp64fc* pSrcDst, int len) {
			return ippsAddC_64fc_I(val, pSrcDst, len);
		}

		inline auto AddC(const Ipp64fc* pSrc, Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsAddC_64fc(pSrc, val, pDst, len);
		}

		inline auto Add(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_8u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsAdd_8u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp16u* pDst, int len) {
			return ippsAdd_8u16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add_16s(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len) {
			return ippsAdd_16s_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len) {
			return ippsAdd_16s(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_16s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsAdd_16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add_16s32s(const Ipp16s* pSrc, Ipp32s* pSrcDst, int len) {
			return ippsAdd_16s32s_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp32f* pDst, int len) {
			return ippsAdd_16s32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add(const Ipp16sc* pSrc, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_16sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp16sc* pSrc1, const Ipp16sc* pSrc2, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsAdd_16sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_16u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsAdd_16u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len) {
			return ippsAdd_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_32s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp32s* pSrc1, const Ipp32s* pSrc2, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsAdd_32s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add(const Ipp32sc* pSrc, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsAdd_32sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Add(const Ipp32sc* pSrc1, const Ipp32sc* pSrc2, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsAdd_32sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add(const Ipp32u* pSrc1, const Ipp32u* pSrc2, Ipp32u* pDst, int len) {
			return ippsAdd_32u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add_32u(const Ipp32u* pSrc, Ipp32u* pSrcDst, int len) {
			return ippsAdd_32u_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp64s* pSrc1, const Ipp64s* pSrc2, Ipp64s* pDst, int len, int scaleFactor) {
			return ippsAdd_64s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Add_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsAdd_32f_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsAdd_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add_32fc(const Ipp32fc* pSrc, Ipp32fc* pSrcDst, int len) {
			return ippsAdd_32fc_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pDst, int len) {
			return ippsAdd_32fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsAdd_64f_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsAdd_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Add_64fc(const Ipp64fc* pSrc, Ipp64fc* pSrcDst, int len) {
			return ippsAdd_64fc_I(pSrc, pSrcDst, len);
		}

		inline auto Add(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, Ipp64fc* pDst, int len) {
			return ippsAdd_64fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto AddProductC(const Ipp32f* pSrc, const Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsAddProductC_32f(pSrc, val, pSrcDst, len);
		}

		inline auto AddProductC(const Ipp64f* pSrc, const Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsAddProductC_64f(pSrc, val, pSrcDst, len);
		}

		inline auto AddProduct(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsAddProduct_16s_Sfs(pSrc1, pSrc2, pSrcDst, len, scaleFactor);
		}

		inline auto AddProduct(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsAddProduct_16s32s_Sfs(pSrc1, pSrc2, pSrcDst, len, scaleFactor);
		}

		inline auto AddProduct(const Ipp32s* pSrc1, const Ipp32s* pSrc2, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsAddProduct_32s_Sfs(pSrc1, pSrc2, pSrcDst, len, scaleFactor);
		}

		inline auto AddProduct(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pSrcDst, int len) {
			return ippsAddProduct_32f(pSrc1, pSrc2, pSrcDst, len);
		}

		inline auto AddProduct(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pSrcDst, int len) {
			return ippsAddProduct_32fc(pSrc1, pSrc2, pSrcDst, len);
		}

		inline auto AddProduct(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pSrcDst, int len) {
			return ippsAddProduct_64f(pSrc1, pSrc2, pSrcDst, len);
		}

		inline auto AddProduct(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, Ipp64fc* pSrcDst, int len) {
			return ippsAddProduct_64fc(pSrc1, pSrc2, pSrcDst, len);
		}

		inline auto MulC(Ipp8u val, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_8u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsMulC_8u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC_16s(Ipp16s val, Ipp16s* pSrcDst, int len) {
			return ippsMulC_16s_I(val, pSrcDst, len);
		}

		inline auto MulC(Ipp16s val, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_16s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp16s* pSrc, Ipp16s val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMulC_16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC(Ipp16sc val, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_16sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp16sc* pSrc, Ipp16sc val, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsMulC_16sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC(Ipp16u val, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_16u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsMulC_16u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC(Ipp32s val, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_32s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp32s* pSrc, Ipp32s val, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsMulC_32s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC(Ipp32sc val, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsMulC_32sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC(const Ipp32sc* pSrc, Ipp32sc val, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsMulC_32sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC(Ipp64s val, Ipp64s* pSrcDst, Ipp32u len, int scaleFactor) {
			return ippsMulC_64s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto MulC_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsMulC_32f_I(val, pSrcDst, len);
		}

		inline auto MulC(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsMulC_32f(pSrc, val, pDst, len);
		}

		inline auto MulC_32fc(Ipp32fc val, Ipp32fc* pSrcDst, int len) {
			return ippsMulC_32fc_I(val, pSrcDst, len);
		}

		inline auto MulC(const Ipp32fc* pSrc, Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsMulC_32fc(pSrc, val, pDst, len);
		}

		inline auto MulC(const Ipp32f* pSrc, Ipp32f val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMulC_32f16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto MulC_Low(const Ipp32f* pSrc, Ipp32f val, Ipp16s* pDst, int len) {
			return ippsMulC_Low_32f16s(pSrc, val, pDst, len);
		}

		inline auto MulC_64f(Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsMulC_64f_I(val, pSrcDst, len);
		}

		inline auto MulC(const Ipp64f* pSrc, Ipp64f val, Ipp64f* pDst, int len) {
			return ippsMulC_64f(pSrc, val, pDst, len);
		}

		inline auto MulC_64fc(Ipp64fc val, Ipp64fc* pSrcDst, int len) {
			return ippsMulC_64fc_I(val, pSrcDst, len);
		}

		inline auto MulC(const Ipp64fc* pSrc, Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsMulC_64fc(pSrc, val, pDst, len);
		}

		inline auto MulC(Ipp64f val, Ipp64s* pSrcDst, Ipp32u len, int scaleFactor) {
			return ippsMulC_64f64s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsMul_8u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsMul_8u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp16u* pDst, int len) {
			return ippsMul_8u16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul_16s(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len) {
			return ippsMul_16s_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len) {
			return ippsMul_16s(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsMul_16s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMul_16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16sc* pSrc, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsMul_16sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16sc* pSrc1, const Ipp16sc* pSrc2, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsMul_16sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsMul_16s32s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp32f* pDst, int len) {
			return ippsMul_16s32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsMul_16u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsMul_16u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp16u* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMul_16u16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsMul_32s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp32s* pSrc1, const Ipp32s* pSrc2, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsMul_32s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp32sc* pSrc, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsMul_32sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Mul(const Ipp32sc* pSrc1, const Ipp32sc* pSrc2, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsMul_32sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Mul_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMul_32f_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsMul_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul_32fc(const Ipp32fc* pSrc, Ipp32fc* pSrcDst, int len) {
			return ippsMul_32fc_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pDst, int len) {
			return ippsMul_32fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul_32f32fc(const Ipp32f* pSrc, Ipp32fc* pSrcDst, int len) {
			return ippsMul_32f32fc_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp32f* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pDst, int len) {
			return ippsMul_32f32fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsMul_64f_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsMul_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Mul_64fc(const Ipp64fc* pSrc, Ipp64fc* pSrcDst, int len) {
			return ippsMul_64fc_I(pSrc, pSrcDst, len);
		}

		inline auto Mul(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, Ipp64fc* pDst, int len) {
			return ippsMul_64fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto SubC(Ipp8u val, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_8u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsSubC_8u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC_16s(Ipp16s val, Ipp16s* pSrcDst, int len) {
			return ippsSubC_16s_I(val, pSrcDst, len);
		}

		inline auto SubC(Ipp16s val, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_16s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp16s* pSrc, Ipp16s val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSubC_16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC(Ipp16sc val, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_16sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp16sc* pSrc, Ipp16sc val, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsSubC_16sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC(Ipp16u val, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_16u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsSubC_16u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC(Ipp32s val, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_32s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp32s* pSrc, Ipp32s val, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsSubC_32s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC(Ipp32sc val, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsSubC_32sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubC(const Ipp32sc* pSrc, Ipp32sc val, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsSubC_32sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubC_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsSubC_32f_I(val, pSrcDst, len);
		}

		inline auto SubC(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsSubC_32f(pSrc, val, pDst, len);
		}

		inline auto SubC_32fc(Ipp32fc val, Ipp32fc* pSrcDst, int len) {
			return ippsSubC_32fc_I(val, pSrcDst, len);
		}

		inline auto SubC(const Ipp32fc* pSrc, Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsSubC_32fc(pSrc, val, pDst, len);
		}

		inline auto SubC_64f(Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsSubC_64f_I(val, pSrcDst, len);
		}

		inline auto SubC(const Ipp64f* pSrc, Ipp64f val, Ipp64f* pDst, int len) {
			return ippsSubC_64f(pSrc, val, pDst, len);
		}

		inline auto SubC_64fc(Ipp64fc val, Ipp64fc* pSrcDst, int len) {
			return ippsSubC_64fc_I(val, pSrcDst, len);
		}

		inline auto SubC(const Ipp64fc* pSrc, Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsSubC_64fc(pSrc, val, pDst, len);
		}

		inline auto SubCRev(Ipp8u val, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_8u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsSubCRev_8u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev(Ipp16s val, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_16s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp16s* pSrc, Ipp16s val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSubCRev_16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev(Ipp16sc val, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_16sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp16sc* pSrc, Ipp16sc val, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsSubCRev_16sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev(Ipp16u val, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_16u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsSubCRev_16u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev(Ipp32s val, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_32s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp32s* pSrc, Ipp32s val, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsSubCRev_32s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev(Ipp32sc val, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsSubCRev_32sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto SubCRev(const Ipp32sc* pSrc, Ipp32sc val, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsSubCRev_32sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto SubCRev_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsSubCRev_32f_I(val, pSrcDst, len);
		}

		inline auto SubCRev(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsSubCRev_32f(pSrc, val, pDst, len);
		}

		inline auto SubCRev_32fc(Ipp32fc val, Ipp32fc* pSrcDst, int len) {
			return ippsSubCRev_32fc_I(val, pSrcDst, len);
		}

		inline auto SubCRev(const Ipp32fc* pSrc, Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsSubCRev_32fc(pSrc, val, pDst, len);
		}

		inline auto SubCRev_64f(Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsSubCRev_64f_I(val, pSrcDst, len);
		}

		inline auto SubCRev(const Ipp64f* pSrc, Ipp64f val, Ipp64f* pDst, int len) {
			return ippsSubCRev_64f(pSrc, val, pDst, len);
		}

		inline auto SubCRev_64fc(Ipp64fc val, Ipp64fc* pSrcDst, int len) {
			return ippsSubCRev_64fc_I(val, pSrcDst, len);
		}

		inline auto SubCRev(const Ipp64fc* pSrc, Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsSubCRev_64fc(pSrc, val, pDst, len);
		}

		inline auto Sub(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsSub_8u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsSub_8u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub_16s(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len) {
			return ippsSub_16s_I(pSrc, pSrcDst, len);
		}

		inline auto Sub(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len) {
			return ippsSub_16s(pSrc1, pSrc2, pDst, len);
		}

		inline auto Sub(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsSub_16s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSub_16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp16sc* pSrc, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsSub_16sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp16sc* pSrc1, const Ipp16sc* pSrc2, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsSub_16sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp32f* pDst, int len) {
			return ippsSub_16s32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Sub(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsSub_16u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsSub_16u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsSub_32s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp32s* pSrc1, const Ipp32s* pSrc2, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsSub_32s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp32sc* pSrc, Ipp32sc* pSrcDst, int len, int scaleFactor) {
			return ippsSub_32sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Sub(const Ipp32sc* pSrc1, const Ipp32sc* pSrc2, Ipp32sc* pDst, int len, int scaleFactor) {
			return ippsSub_32sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Sub_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsSub_32f_I(pSrc, pSrcDst, len);
		}

		inline auto Sub(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsSub_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Sub_32fc(const Ipp32fc* pSrc, Ipp32fc* pSrcDst, int len) {
			return ippsSub_32fc_I(pSrc, pSrcDst, len);
		}

		inline auto Sub(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pDst, int len) {
			return ippsSub_32fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Sub_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsSub_64f_I(pSrc, pSrcDst, len);
		}

		inline auto Sub(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsSub_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Sub_64fc(const Ipp64fc* pSrc, Ipp64fc* pSrcDst, int len) {
			return ippsSub_64fc_I(pSrc, pSrcDst, len);
		}

		inline auto Sub(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, Ipp64fc* pDst, int len) {
			return ippsSub_64fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto DivC(Ipp8u val, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsDivC_8u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto DivC(const Ipp8u* pSrc, Ipp8u val, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsDivC_8u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto DivC(Ipp16s val, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsDivC_16s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto DivC(const Ipp16s* pSrc, Ipp16s val, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsDivC_16s_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto DivC(Ipp16sc val, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsDivC_16sc_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto DivC(const Ipp16sc* pSrc, Ipp16sc val, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsDivC_16sc_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto DivC(Ipp16u val, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsDivC_16u_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto DivC(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsDivC_16u_Sfs(pSrc, val, pDst, len, scaleFactor);
		}

		inline auto DivC(Ipp64s val, Ipp64s* pSrcDst, Ipp32u len, int scaleFactor) {
			return ippsDivC_64s_ISfs(val, pSrcDst, len, scaleFactor);
		}

		inline auto DivC_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsDivC_32f_I(val, pSrcDst, len);
		}

		inline auto DivC(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsDivC_32f(pSrc, val, pDst, len);
		}

		inline auto DivC_32fc(Ipp32fc val, Ipp32fc* pSrcDst, int len) {
			return ippsDivC_32fc_I(val, pSrcDst, len);
		}

		inline auto DivC(const Ipp32fc* pSrc, Ipp32fc val, Ipp32fc* pDst, int len) {
			return ippsDivC_32fc(pSrc, val, pDst, len);
		}

		inline auto DivC_64f(Ipp64f val, Ipp64f* pSrcDst, int len) {
			return ippsDivC_64f_I(val, pSrcDst, len);
		}

		inline auto DivC(const Ipp64f* pSrc, Ipp64f val, Ipp64f* pDst, int len) {
			return ippsDivC_64f(pSrc, val, pDst, len);
		}

		inline auto DivC_64fc(Ipp64fc val, Ipp64fc* pSrcDst, int len) {
			return ippsDivC_64fc_I(val, pSrcDst, len);
		}

		inline auto DivC(const Ipp64fc* pSrc, Ipp64fc val, Ipp64fc* pDst, int len) {
			return ippsDivC_64fc(pSrc, val, pDst, len);
		}

		inline auto DivCRev_16u(Ipp16u val, Ipp16u* pSrcDst, int len) {
			return ippsDivCRev_16u_I(val, pSrcDst, len);
		}

		inline auto DivCRev(const Ipp16u* pSrc, Ipp16u val, Ipp16u* pDst, int len) {
			return ippsDivCRev_16u(pSrc, val, pDst, len);
		}

		inline auto DivCRev_32f(Ipp32f val, Ipp32f* pSrcDst, int len) {
			return ippsDivCRev_32f_I(val, pSrcDst, len);
		}

		inline auto DivCRev(const Ipp32f* pSrc, Ipp32f val, Ipp32f* pDst, int len) {
			return ippsDivCRev_32f(pSrc, val, pDst, len);
		}

		inline auto Div(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsDiv_8u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Div(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsDiv_8u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsDiv_16s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsDiv_16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16sc* pSrc, Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsDiv_16sc_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16sc* pSrc1, const Ipp16sc* pSrc2, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsDiv_16sc_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsDiv_32s_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Div(const Ipp32s* pSrc1, const Ipp32s* pSrc2, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsDiv_32s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16s* pSrc1, const Ipp32s* pSrc2, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsDiv_32s16s_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsDiv_16u_ISfs(pSrc, pSrcDst, len, scaleFactor);
		}

		inline auto Div(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsDiv_16u_Sfs(pSrc1, pSrc2, pDst, len, scaleFactor);
		}

		inline auto Div_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsDiv_32f_I(pSrc, pSrcDst, len);
		}

		inline auto Div(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsDiv_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Div_32fc(const Ipp32fc* pSrc, Ipp32fc* pSrcDst, int len) {
			return ippsDiv_32fc_I(pSrc, pSrcDst, len);
		}

		inline auto Div(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pDst, int len) {
			return ippsDiv_32fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Div_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsDiv_64f_I(pSrc, pSrcDst, len);
		}

		inline auto Div(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsDiv_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto Div_64fc(const Ipp64fc* pSrc, Ipp64fc* pSrcDst, int len) {
			return ippsDiv_64fc_I(pSrc, pSrcDst, len);
		}

		inline auto Div(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, Ipp64fc* pDst, int len) {
			return ippsDiv_64fc(pSrc1, pSrc2, pDst, len);
		}

		inline auto Div_Round(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_8u_ISfs(pSrc, pSrcDst, len, rndMode, scaleFactor);
		}

		inline auto Div_Round(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_8u_Sfs(pSrc1, pSrc2, pDst, len, rndMode, scaleFactor);
		}

		inline auto Div_Round(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_16s_ISfs(pSrc, pSrcDst, len, rndMode, scaleFactor);
		}

		inline auto Div_Round(const Ipp16s* pSrc1, const Ipp16s* pSrc2, Ipp16s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_16s_Sfs(pSrc1, pSrc2, pDst, len, rndMode, scaleFactor);
		}

		inline auto Div_Round(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_16u_ISfs(pSrc, pSrcDst, len, rndMode, scaleFactor);
		}

		inline auto Div_Round(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsDiv_Round_16u_Sfs(pSrc1, pSrc2, pDst, len, rndMode, scaleFactor);
		}

		inline auto Abs_16s(Ipp16s* pSrcDst, int len) {
			return ippsAbs_16s_I(pSrcDst, len);
		}

		inline auto Abs(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsAbs_16s(pSrc, pDst, len);
		}

		inline auto Abs_32s(Ipp32s* pSrcDst, int len) {
			return ippsAbs_32s_I(pSrcDst, len);
		}

		inline auto Abs(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
			return ippsAbs_32s(pSrc, pDst, len);
		}

		inline auto Abs_32f(Ipp32f* pSrcDst, int len) {
			return ippsAbs_32f_I(pSrcDst, len);
		}

		inline auto Abs(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsAbs_32f(pSrc, pDst, len);
		}

		inline auto Abs_64f(Ipp64f* pSrcDst, int len) {
			return ippsAbs_64f_I(pSrcDst, len);
		}

		inline auto Abs(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsAbs_64f(pSrc, pDst, len);
		}

		inline auto Sqr(Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsSqr_8u_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqr(const Ipp8u* pSrc, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsSqr_8u_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqr(Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsSqr_16s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqr(const Ipp16s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSqr_16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqr(Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsSqr_16sc_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqr(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsSqr_16sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqr(Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsSqr_16u_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqr(const Ipp16u* pSrc, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsSqr_16u_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqr_32f(Ipp32f* pSrcDst, int len) {
			return ippsSqr_32f_I(pSrcDst, len);
		}

		inline auto Sqr(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsSqr_32f(pSrc, pDst, len);
		}

		inline auto Sqr_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsSqr_32fc_I(pSrcDst, len);
		}

		inline auto Sqr(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsSqr_32fc(pSrc, pDst, len);
		}

		inline auto Sqr_64f(Ipp64f* pSrcDst, int len) {
			return ippsSqr_64f_I(pSrcDst, len);
		}

		inline auto Sqr(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsSqr_64f(pSrc, pDst, len);
		}

		inline auto Sqr_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsSqr_64fc_I(pSrcDst, len);
		}

		inline auto Sqr(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsSqr_64fc(pSrc, pDst, len);
		}

		inline auto Sqrt(Ipp8u* pSrcDst, int len, int scaleFactor) {
			return ippsSqrt_8u_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqrt(const Ipp8u* pSrc, Ipp8u* pDst, int len, int scaleFactor) {
			return ippsSqrt_8u_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqrt(Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsSqrt_16s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqrt(const Ipp16s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSqrt_16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqrt(Ipp16sc* pSrcDst, int len, int scaleFactor) {
			return ippsSqrt_16sc_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqrt(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, int scaleFactor) {
			return ippsSqrt_16sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqrt(Ipp16u* pSrcDst, int len, int scaleFactor) {
			return ippsSqrt_16u_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Sqrt(const Ipp16u* pSrc, Ipp16u* pDst, int len, int scaleFactor) {
			return ippsSqrt_16u_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqrt(const Ipp32s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsSqrt_32s16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Sqrt_32f(Ipp32f* pSrcDst, int len) {
			return ippsSqrt_32f_I(pSrcDst, len);
		}

		inline auto Sqrt(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsSqrt_32f(pSrc, pDst, len);
		}

		inline auto Sqrt_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsSqrt_32fc_I(pSrcDst, len);
		}

		inline auto Sqrt(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsSqrt_32fc(pSrc, pDst, len);
		}

		inline auto Sqrt_64f(Ipp64f* pSrcDst, int len) {
			return ippsSqrt_64f_I(pSrcDst, len);
		}

		inline auto Sqrt(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsSqrt_64f(pSrc, pDst, len);
		}

		inline auto Sqrt_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsSqrt_64fc_I(pSrcDst, len);
		}

		inline auto Sqrt(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsSqrt_64fc(pSrc, pDst, len);
		}

		inline auto Cubrt(const Ipp32s* pSrc, Ipp16s* pDst, int Len, int scaleFactor) {
			return ippsCubrt_32s16s_Sfs(pSrc, pDst, Len, scaleFactor);
		}

		inline auto Cubrt(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsCubrt_32f(pSrc, pDst, len);
		}

		inline auto Exp(Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsExp_16s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Exp(const Ipp16s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsExp_16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Exp(Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsExp_32s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Exp(const Ipp32s* pSrc, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsExp_32s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Exp_32f(Ipp32f* pSrcDst, int len) {
			return ippsExp_32f_I(pSrcDst, len);
		}

		inline auto Exp(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsExp_32f(pSrc, pDst, len);
		}

		inline auto Exp_64f(Ipp64f* pSrcDst, int len) {
			return ippsExp_64f_I(pSrcDst, len);
		}

		inline auto Exp(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsExp_64f(pSrc, pDst, len);
		}

		inline auto Ln(Ipp16s* pSrcDst, int len, int scaleFactor) {
			return ippsLn_16s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Ln(const Ipp16s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsLn_16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Ln(Ipp32s* pSrcDst, int len, int scaleFactor) {
			return ippsLn_32s_ISfs(pSrcDst, len, scaleFactor);
		}

		inline auto Ln(const Ipp32s* pSrc, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsLn_32s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Ln_32f(Ipp32f* pSrcDst, int len) {
			return ippsLn_32f_I(pSrcDst, len);
		}

		inline auto Ln(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsLn_32f(pSrc, pDst, len);
		}

		inline auto Ln_64f(Ipp64f* pSrcDst, int len) {
			return ippsLn_64f_I(pSrcDst, len);
		}

		inline auto Ln(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsLn_64f(pSrc, pDst, len);
		}

		inline auto SumLn(const Ipp16s* pSrc, int len, Ipp32f* pSum) {
			return ippsSumLn_16s32f(pSrc, len, pSum);
		}

		inline auto SumLn(const Ipp32f* pSrc, int len, Ipp32f* pSum) {
			return ippsSumLn_32f(pSrc, len, pSum);
		}

		inline auto SumLn(const Ipp32f* pSrc, int len, Ipp64f* pSum) {
			return ippsSumLn_32f64f(pSrc, len, pSum);
		}

		inline auto SumLn(const Ipp64f* pSrc, int len, Ipp64f* pSum) {
			return ippsSumLn_64f(pSrc, len, pSum);
		}

		inline auto Arctan_32f(Ipp32f* pSrcDst, int len) {
			return ippsArctan_32f_I(pSrcDst, len);
		}

		inline auto Arctan(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsArctan_32f(pSrc, pDst, len);
		}

		inline auto Arctan_64f(Ipp64f* pSrcDst, int len) {
			return ippsArctan_64f_I(pSrcDst, len);
		}

		inline auto Arctan(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsArctan_64f(pSrc, pDst, len);
		}

		inline auto Normalize(Ipp16s* pSrcDst, int len, Ipp16s vSub, int vDiv, int scaleFactor) {
			return ippsNormalize_16s_ISfs(pSrcDst, len, vSub, vDiv, scaleFactor);
		}

		inline auto Normalize(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s vSub, int vDiv, int scaleFactor) {
			return ippsNormalize_16s_Sfs(pSrc, pDst, len, vSub, vDiv, scaleFactor);
		}

		inline auto Normalize(Ipp16sc* pSrcDst, int len, Ipp16sc vSub, int vDiv, int scaleFactor) {
			return ippsNormalize_16sc_ISfs(pSrcDst, len, vSub, vDiv, scaleFactor);
		}

		inline auto Normalize(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16sc vSub, int vDiv, int scaleFactor) {
			return ippsNormalize_16sc_Sfs(pSrc, pDst, len, vSub, vDiv, scaleFactor);
		}

		inline auto Normalize_32f(Ipp32f* pSrcDst, int len, Ipp32f vSub, Ipp32f vDiv) {
			return ippsNormalize_32f_I(pSrcDst, len, vSub, vDiv);
		}

		inline auto Normalize(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f vSub, Ipp32f vDiv) {
			return ippsNormalize_32f(pSrc, pDst, len, vSub, vDiv);
		}

		inline auto Normalize_32fc(Ipp32fc* pSrcDst, int len, Ipp32fc vSub, Ipp32f vDiv) {
			return ippsNormalize_32fc_I(pSrcDst, len, vSub, vDiv);
		}

		inline auto Normalize(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32fc vSub, Ipp32f vDiv) {
			return ippsNormalize_32fc(pSrc, pDst, len, vSub, vDiv);
		}

		inline auto Normalize_64f(Ipp64f* pSrcDst, int len, Ipp64f vSub, Ipp64f vDiv) {
			return ippsNormalize_64f_I(pSrcDst, len, vSub, vDiv);
		}

		inline auto Normalize(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f vSub, Ipp64f vDiv) {
			return ippsNormalize_64f(pSrc, pDst, len, vSub, vDiv);
		}

		inline auto Normalize_64fc(Ipp64fc* pSrcDst, int len, Ipp64fc vSub, Ipp64f vDiv) {
			return ippsNormalize_64fc_I(pSrcDst, len, vSub, vDiv);
		}

		inline auto Normalize(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64fc vSub, Ipp64f vDiv) {
			return ippsNormalize_64fc(pSrc, pDst, len, vSub, vDiv);
		}

		inline auto SortAscend_8u(Ipp8u* pSrcDst, int len) {
			return ippsSortAscend_8u_I(pSrcDst, len);
		}

		inline auto SortAscend_16s(Ipp16s* pSrcDst, int len) {
			return ippsSortAscend_16s_I(pSrcDst, len);
		}

		inline auto SortAscend_16u(Ipp16u* pSrcDst, int len) {
			return ippsSortAscend_16u_I(pSrcDst, len);
		}

		inline auto SortAscend_32s(Ipp32s* pSrcDst, int len) {
			return ippsSortAscend_32s_I(pSrcDst, len);
		}

		inline auto SortAscend_32f(Ipp32f* pSrcDst, int len) {
			return ippsSortAscend_32f_I(pSrcDst, len);
		}

		inline auto SortAscend_64f(Ipp64f* pSrcDst, int len) {
			return ippsSortAscend_64f_I(pSrcDst, len);
		}

		inline auto SortDescend_8u(Ipp8u* pSrcDst, int len) {
			return ippsSortDescend_8u_I(pSrcDst, len);
		}

		inline auto SortDescend_16s(Ipp16s* pSrcDst, int len) {
			return ippsSortDescend_16s_I(pSrcDst, len);
		}

		inline auto SortDescend_16u(Ipp16u* pSrcDst, int len) {
			return ippsSortDescend_16u_I(pSrcDst, len);
		}

		inline auto SortDescend_32s(Ipp32s* pSrcDst, int len) {
			return ippsSortDescend_32s_I(pSrcDst, len);
		}

		inline auto SortDescend_32f(Ipp32f* pSrcDst, int len) {
			return ippsSortDescend_32f_I(pSrcDst, len);
		}

		inline auto SortDescend_64f(Ipp64f* pSrcDst, int len) {
			return ippsSortDescend_64f_I(pSrcDst, len);
		}

		inline auto SortIndexAscend_8u(Ipp8u* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_8u_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexAscend_16s(Ipp16s* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_16s_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexAscend_16u(Ipp16u* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_16u_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexAscend_32s(Ipp32s* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_32s_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexAscend_32f(Ipp32f* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_32f_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexAscend_64f(Ipp64f* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexAscend_64f_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_8u(Ipp8u* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_8u_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_16s(Ipp16s* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_16s_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_16u(Ipp16u* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_16u_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_32s(Ipp32s* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_32s_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_32f(Ipp32f* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_32f_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortIndexDescend_64f(Ipp64f* pSrcDst, int* pDstIdx, int len) {
			return ippsSortIndexDescend_64f_I(pSrcDst, pDstIdx, len);
		}

		inline auto SortRadixAscend_8u(Ipp8u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_8u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_16u(Ipp16u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_16u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_16s(Ipp16s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_16s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_32u(Ipp32u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_32u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_32s(Ipp32s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_32s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_32f(Ipp32f* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_32f_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64u(Ipp64u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64s(Ipp64s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64f(Ipp64f* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64f_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_8u(Ipp8u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_8u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_16u(Ipp16u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_16u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_16s(Ipp16s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_16s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_32u(Ipp32u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_32u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_32s(Ipp32s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_32s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_32f(Ipp32f* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_32f_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64u(Ipp64u* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64u_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64s(Ipp64s* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64s_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64f(Ipp64f* pSrcDst, int len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64f_I(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp8u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_8u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp16s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_16s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp16u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_16u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp32s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_32s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp32u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_32u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp32f* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_32f(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp64s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_64s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp64u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_64u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend(const Ipp64f* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_64f(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp8u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_8u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp16s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_16s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp16u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_16u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp32s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_32s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp32u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_32u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp32f* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_32f(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp64s* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_64s(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp64u* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_64u(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend(const Ipp64f* pSrc, Ipp32s srcStrideBytes, Ipp32s* pDstIndx, int len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_64f(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SwapBytes_16u(Ipp16u* pSrcDst, int len) {
			return ippsSwapBytes_16u_I(pSrcDst, len);
		}

		inline auto SwapBytes(const Ipp16u* pSrc, Ipp16u* pDst, int len) {
			return ippsSwapBytes_16u(pSrc, pDst, len);
		}

		inline auto SwapBytes_24u(Ipp24u* pSrcDst, int len) {
			return ippsSwapBytes_24u_I(reinterpret_cast<Ipp8u*>(pSrcDst), len);
		}

		inline auto SwapBytes(const Ipp24u* pSrc, Ipp24u* pDst, int len) {
			return ippsSwapBytes_24u(reinterpret_cast<const Ipp8u*>(pSrc), reinterpret_cast<Ipp8u*>(pDst), len);
		}

		inline auto SwapBytes_32u(Ipp32u* pSrcDst, int len) {
			return ippsSwapBytes_32u_I(pSrcDst, len);
		}

		inline auto SwapBytes(const Ipp32u* pSrc, Ipp32u* pDst, int len) {
			return ippsSwapBytes_32u(pSrc, pDst, len);
		}

		inline auto SwapBytes_64u(Ipp64u* pSrcDst, int len) {
			return ippsSwapBytes_64u_I(pSrcDst, len);
		}

		inline auto SwapBytes(const Ipp64u* pSrc, Ipp64u* pDst, int len) {
			return ippsSwapBytes_64u(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp8s* pSrc, Ipp16s* pDst, int len) {
			return ippsConvert_8s16s(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp8s* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_8s32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp8u* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_8u32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp8s* pDst, Ipp32u len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_16s8s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp32s* pDst, int len) {
			return ippsConvert_16s32s(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp16f* pDst, int len, IppRoundMode rndMode) {
			return ippsConvert_16s16f(pSrc, reinterpret_cast<::Ipp16f*>(pDst), len, rndMode);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_16s32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp32f* pDst, int len, int scaleFactor) {
			return ippsConvert_16s32f_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Convert(const Ipp16s* pSrc, Ipp64f* pDst, int len, int scaleFactor) {
			return ippsConvert_16s64f_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Convert(const Ipp16u* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_16u32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp24u* pSrc, Ipp32s* pDst, int len) {
			return ippsConvert_24s32s(reinterpret_cast<const Ipp8u*>(pSrc), pDst, len);
		}

		inline auto Convert(const Ipp24s* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_24s32f(reinterpret_cast<const Ipp8u*>(pSrc), pDst, len);
		}

		inline auto Convert(const Ipp24u* pSrc, Ipp32u* pDst, int len) {
			return ippsConvert_24u32u(reinterpret_cast<const Ipp8u*>(pSrc), pDst, len);
		}

		inline auto Convert(const Ipp24u* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_24u32f(reinterpret_cast<const Ipp8u*>(pSrc), pDst, len);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp16s* pDst, int len) {
			return ippsConvert_32s16s(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsConvert_32s16s_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp24s* pDst, int len, int scaleFactor) {
			return ippsConvert_32s24s_Sfs(pSrc, reinterpret_cast<Ipp8u*>(pDst), len, scaleFactor);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_32s32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp32f* pDst, int len, int scaleFactor) {
			return ippsConvert_32s32f_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp64f* pDst, int len) {
			return ippsConvert_32s64f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp32s* pSrc, Ipp64f* pDst, int len, int scaleFactor) {
			return ippsConvert_32s64f_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Convert(const Ipp32u* pSrc, Ipp24u* pDst, int len, int scaleFactor) {
			return ippsConvert_32u24u_Sfs(pSrc, reinterpret_cast<Ipp8u*>(pDst), len, scaleFactor);
		}

		inline auto Convert(const Ipp64s* pSrc, Ipp32s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64s32s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64s* pSrc, Ipp64f* pDst, Ipp32u len) {
			return ippsConvert_64s64f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp16f* pSrc, Ipp16s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_16f16s_Sfs(reinterpret_cast<const ::Ipp16f*>(pDst), pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp16f* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_16f32f(reinterpret_cast<const ::Ipp16f*>(pDst), pDst, len);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp8s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_32f8s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp8u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_32f8u_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp16s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_32f16s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp16u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_32f16u_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp24s* pDst, int len, int scaleFactor) {
			return ippsConvert_32f24s_Sfs(pSrc, reinterpret_cast<Ipp8u*>(pDst), len, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp24u* pDst, int len, int scaleFactor) {
			return ippsConvert_32f24u_Sfs(pSrc, reinterpret_cast<Ipp8u*>(pDst), len, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp32s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_32f32s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp16f* pDst, int len, IppRoundMode rndMode) {
			return ippsConvert_32f16f(pSrc, reinterpret_cast<::Ipp16f*>(pDst), len, rndMode);
		}

		inline auto Convert(const Ipp32f* pSrc, Ipp64f* pDst, int len) {
			return ippsConvert_32f64f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp16s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f16s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp32s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f32s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp64s* pDst, Ipp32u len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f64s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp32f* pDst, int len) {
			return ippsConvert_64f32f(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp8s* pSrc, Ipp8u* pDst, int len) {
			return ippsConvert_8s8u(pSrc, pDst, len);
		}

		inline auto Convert(const Ipp8u* pSrc, Ipp8s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_8u8s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp8s* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f8s_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp8u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f8u_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Convert(const Ipp64f* pSrc, Ipp16u* pDst, int len, IppRoundMode rndMode, int scaleFactor) {
			return ippsConvert_64f16u_Sfs(pSrc, pDst, len, rndMode, scaleFactor);
		}

		inline auto Conj_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsConj_16sc_I(pSrcDst, len);
		}

		inline auto Conj(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsConj_16sc(pSrc, pDst, len);
		}

		inline auto Conj_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsConj_32fc_I(pSrcDst, len);
		}

		inline auto Conj(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsConj_32fc(pSrc, pDst, len);
		}

		inline auto Conj_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsConj_64fc_I(pSrcDst, len);
		}

		inline auto Conj(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsConj_64fc(pSrc, pDst, len);
		}

		inline auto ConjFlip(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsConjFlip_16sc(pSrc, pDst, len);
		}

		inline auto ConjFlip(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsConjFlip_32fc(pSrc, pDst, len);
		}

		inline auto ConjFlip(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsConjFlip_64fc(pSrc, pDst, len);
		}

		inline auto ConjCcs_32fc(Ipp32fc* pSrcDst, int lenDst) {
			return ippsConjCcs_32fc_I(pSrcDst, lenDst);
		}

		inline auto ConjCcs(const Ipp32f* pSrc, Ipp32fc* pDst, int lenDst) {
			return ippsConjCcs_32fc(pSrc, pDst, lenDst);
		}

		inline auto ConjCcs_64fc(Ipp64fc* pSrcDst, int lenDst) {
			return ippsConjCcs_64fc_I(pSrcDst, lenDst);
		}

		inline auto ConjCcs(const Ipp64f* pSrc, Ipp64fc* pDst, int lenDst) {
			return ippsConjCcs_64fc(pSrc, pDst, lenDst);
		}

		inline auto ConjPack_32fc(Ipp32fc* pSrcDst, int lenDst) {
			return ippsConjPack_32fc_I(pSrcDst, lenDst);
		}

		inline auto ConjPack(const Ipp32f* pSrc, Ipp32fc* pDst, int lenDst) {
			return ippsConjPack_32fc(pSrc, pDst, lenDst);
		}

		inline auto ConjPack_64fc(Ipp64fc* pSrcDst, int lenDst) {
			return ippsConjPack_64fc_I(pSrcDst, lenDst);
		}

		inline auto ConjPack(const Ipp64f* pSrc, Ipp64fc* pDst, int lenDst) {
			return ippsConjPack_64fc(pSrc, pDst, lenDst);
		}

		inline auto ConjPerm_32fc(Ipp32fc* pSrcDst, int lenDst) {
			return ippsConjPerm_32fc_I(pSrcDst, lenDst);
		}

		inline auto ConjPerm(const Ipp32f* pSrc, Ipp32fc* pDst, int lenDst) {
			return ippsConjPerm_32fc(pSrc, pDst, lenDst);
		}

		inline auto ConjPerm_64fc(Ipp64fc* pSrcDst, int lenDst) {
			return ippsConjPerm_64fc_I(pSrcDst, lenDst);
		}

		inline auto ConjPerm(const Ipp64f* pSrc, Ipp64fc* pDst, int lenDst) {
			return ippsConjPerm_64fc(pSrc, pDst, lenDst);
		}

		inline auto Magnitude(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMagnitude_16s_Sfs(pSrcRe, pSrcIm, pDst, len, scaleFactor);
		}

		inline auto Magnitude(const Ipp16sc* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsMagnitude_16sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Magnitude(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp32f* pDst, int len) {
			return ippsMagnitude_16s32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto Magnitude(const Ipp16sc* pSrc, Ipp32f* pDst, int len) {
			return ippsMagnitude_16sc32f(pSrc, pDst, len);
		}

		inline auto Magnitude(const Ipp32sc* pSrc, Ipp32s* pDst, int len, int scaleFactor) {
			return ippsMagnitude_32sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Magnitude(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDst, int len) {
			return ippsMagnitude_32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto Magnitude(const Ipp32fc* pSrc, Ipp32f* pDst, int len) {
			return ippsMagnitude_32fc(pSrc, pDst, len);
		}

		inline auto Magnitude(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDst, int len) {
			return ippsMagnitude_64f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto Magnitude(const Ipp64fc* pSrc, Ipp64f* pDst, int len) {
			return ippsMagnitude_64fc(pSrc, pDst, len);
		}

		inline auto Phase(const Ipp16sc* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsPhase_16sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto Phase(const Ipp16sc* pSrc, Ipp32f* pDst, int len) {
			return ippsPhase_16sc32f(pSrc, pDst, len);
		}

		inline auto Phase(const Ipp64fc* pSrc, Ipp64f* pDst, int len) {
			return ippsPhase_64fc(pSrc, pDst, len);
		}

		inline auto Phase(const Ipp32fc* pSrc, Ipp32f* pDst, int len) {
			return ippsPhase_32fc(pSrc, pDst, len);
		}

		inline auto Phase(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsPhase_16s_Sfs(pSrcRe, pSrcIm, pDst, len, scaleFactor);
		}

		inline auto Phase(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp32f* pDst, int len) {
			return ippsPhase_16s32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto Phase(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDst, int len) {
			return ippsPhase_64f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto Phase(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDst, int len) {
			return ippsPhase_32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto PowerSpectr(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsPowerSpectr_16s_Sfs(pSrcRe, pSrcIm, pDst, len, scaleFactor);
		}

		inline auto PowerSpectr(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp32f* pDst, int len) {
			return ippsPowerSpectr_16s32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto PowerSpectr(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDst, int len) {
			return ippsPowerSpectr_32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto PowerSpectr(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDst, int len) {
			return ippsPowerSpectr_64f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto PowerSpectr(const Ipp16sc* pSrc, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsPowerSpectr_16sc_Sfs(pSrc, pDst, len, scaleFactor);
		}

		inline auto PowerSpectr(const Ipp16sc* pSrc, Ipp32f* pDst, int len) {
			return ippsPowerSpectr_16sc32f(pSrc, pDst, len);
		}

		inline auto PowerSpectr(const Ipp32fc* pSrc, Ipp32f* pDst, int len) {
			return ippsPowerSpectr_32fc(pSrc, pDst, len);
		}

		inline auto PowerSpectr(const Ipp64fc* pSrc, Ipp64f* pDst, int len) {
			return ippsPowerSpectr_64fc(pSrc, pDst, len);
		}

		inline auto Real(const Ipp64fc* pSrc, Ipp64f* pDstRe, int len) {
			return ippsReal_64fc(pSrc, pDstRe, len);
		}

		inline auto Real(const Ipp32fc* pSrc, Ipp32f* pDstRe, int len) {
			return ippsReal_32fc(pSrc, pDstRe, len);
		}

		inline auto Real(const Ipp16sc* pSrc, Ipp16s* pDstRe, int len) {
			return ippsReal_16sc(pSrc, pDstRe, len);
		}

		inline auto Imag(const Ipp64fc* pSrc, Ipp64f* pDstIm, int len) {
			return ippsImag_64fc(pSrc, pDstIm, len);
		}

		inline auto Imag(const Ipp32fc* pSrc, Ipp32f* pDstIm, int len) {
			return ippsImag_32fc(pSrc, pDstIm, len);
		}

		inline auto Imag(const Ipp16sc* pSrc, Ipp16s* pDstIm, int len) {
			return ippsImag_16sc(pSrc, pDstIm, len);
		}

		inline auto RealToCplx(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64fc* pDst, int len) {
			return ippsRealToCplx_64f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto RealToCplx(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32fc* pDst, int len) {
			return ippsRealToCplx_32f(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto RealToCplx(const Ipp16s* pSrcRe, const Ipp16s* pSrcIm, Ipp16sc* pDst, int len) {
			return ippsRealToCplx_16s(pSrcRe, pSrcIm, pDst, len);
		}

		inline auto CplxToReal(const Ipp64fc* pSrc, Ipp64f* pDstRe, Ipp64f* pDstIm, int len) {
			return ippsCplxToReal_64fc(pSrc, pDstRe, pDstIm, len);
		}

		inline auto CplxToReal(const Ipp32fc* pSrc, Ipp32f* pDstRe, Ipp32f* pDstIm, int len) {
			return ippsCplxToReal_32fc(pSrc, pDstRe, pDstIm, len);
		}

		inline auto CplxToReal(const Ipp16sc* pSrc, Ipp16s* pDstRe, Ipp16s* pDstIm, int len) {
			return ippsCplxToReal_16sc(pSrc, pDstRe, pDstIm, len);
		}

		inline auto Threshold_16s(Ipp16s* pSrcDst, int len, Ipp16s level, IppCmpOp relOp) {
			return ippsThreshold_16s_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold_16sc(Ipp16sc* pSrcDst, int len, Ipp16s level, IppCmpOp relOp) {
			return ippsThreshold_16sc_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold_32f(Ipp32f* pSrcDst, int len, Ipp32f level, IppCmpOp relOp) {
			return ippsThreshold_32f_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level, IppCmpOp relOp) {
			return ippsThreshold_32fc_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold_64f(Ipp64f* pSrcDst, int len, Ipp64f level, IppCmpOp relOp) {
			return ippsThreshold_64f_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level, IppCmpOp relOp) {
			return ippsThreshold_64fc_I(pSrcDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level, IppCmpOp relOp) {
			return ippsThreshold_16s(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16s level, IppCmpOp relOp) {
			return ippsThreshold_16sc(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level, IppCmpOp relOp) {
			return ippsThreshold_32f(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level, IppCmpOp relOp) {
			return ippsThreshold_32fc(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level, IppCmpOp relOp) {
			return ippsThreshold_64f(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level, IppCmpOp relOp) {
			return ippsThreshold_64fc(pSrc, pDst, len, level, relOp);
		}

		inline auto Threshold_LT_16s(Ipp16s* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_LT_16s_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_16sc(Ipp16sc* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_LT_16sc_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_32s(Ipp32s* pSrcDst, int len, Ipp32s level) {
			return ippsThreshold_LT_32s_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_32f(Ipp32f* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_LT_32f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_LT_32fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_64f(Ipp64f* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_LT_64f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_LT_64fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_LT(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level) {
			return ippsThreshold_LT_16s(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16s level) {
			return ippsThreshold_LT_16sc(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s level) {
			return ippsThreshold_LT_32s(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level) {
			return ippsThreshold_LT_32f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level) {
			return ippsThreshold_LT_32fc(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level) {
			return ippsThreshold_LT_64f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LT(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level) {
			return ippsThreshold_LT_64fc(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT_16s(Ipp16s* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_GT_16s_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_16sc(Ipp16sc* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_GT_16sc_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_32s(Ipp32s* pSrcDst, int len, Ipp32s level) {
			return ippsThreshold_GT_32s_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_32f(Ipp32f* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_GT_32f_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_GT_32fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_64f(Ipp64f* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_GT_64f_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_GT_64fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_GT(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level) {
			return ippsThreshold_GT_16s(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16s level) {
			return ippsThreshold_GT_16sc(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s level) {
			return ippsThreshold_GT_32s(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level) {
			return ippsThreshold_GT_32f(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level) {
			return ippsThreshold_GT_32fc(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level) {
			return ippsThreshold_GT_64f(pSrc, pDst, len, level);
		}

		inline auto Threshold_GT(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level) {
			return ippsThreshold_GT_64fc(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTAbs(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level) {
			return ippsThreshold_LTAbs_16s(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTAbs(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s level) {
			return ippsThreshold_LTAbs_32s(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTAbs(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level) {
			return ippsThreshold_LTAbs_32f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTAbs(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level) {
			return ippsThreshold_LTAbs_64f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTAbs_16s(Ipp16s* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_LTAbs_16s_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTAbs_32s(Ipp32s* pSrcDst, int len, Ipp32s level) {
			return ippsThreshold_LTAbs_32s_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTAbs_32f(Ipp32f* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_LTAbs_32f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTAbs_64f(Ipp64f* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_LTAbs_64f_I(pSrcDst, len, level);
		}

		inline auto Threshold_GTAbs(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level) {
			return ippsThreshold_GTAbs_16s(pSrc, pDst, len, level);
		}

		inline auto Threshold_GTAbs(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s level) {
			return ippsThreshold_GTAbs_32s(pSrc, pDst, len, level);
		}

		inline auto Threshold_GTAbs(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level) {
			return ippsThreshold_GTAbs_32f(pSrc, pDst, len, level);
		}

		inline auto Threshold_GTAbs(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level) {
			return ippsThreshold_GTAbs_64f(pSrc, pDst, len, level);
		}

		inline auto Threshold_GTAbs_16s(Ipp16s* pSrcDst, int len, Ipp16s level) {
			return ippsThreshold_GTAbs_16s_I(pSrcDst, len, level);
		}

		inline auto Threshold_GTAbs_32s(Ipp32s* pSrcDst, int len, Ipp32s level) {
			return ippsThreshold_GTAbs_32s_I(pSrcDst, len, level);
		}

		inline auto Threshold_GTAbs_32f(Ipp32f* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_GTAbs_32f_I(pSrcDst, len, level);
		}

		inline auto Threshold_GTAbs_64f(Ipp64f* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_GTAbs_64f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTAbsVal(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_LTAbsVal_16s(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s level, Ipp32s value) {
			return ippsThreshold_LTAbsVal_32s(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_LTAbsVal_32f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_LTAbsVal_64f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal_16s(Ipp16s* pSrcDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_LTAbsVal_16s_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal_32s(Ipp32s* pSrcDst, int len, Ipp32s level, Ipp32s value) {
			return ippsThreshold_LTAbsVal_32s_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal_32f(Ipp32f* pSrcDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_LTAbsVal_32f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTAbsVal_64f(Ipp64f* pSrcDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_LTAbsVal_64f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_16s(Ipp16s* pSrcDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_LTVal_16s_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_16sc(Ipp16sc* pSrcDst, int len, Ipp16s level, Ipp16sc value) {
			return ippsThreshold_LTVal_16sc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_32f(Ipp32f* pSrcDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_LTVal_32f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level, Ipp32fc value) {
			return ippsThreshold_LTVal_32fc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_64f(Ipp64f* pSrcDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_LTVal_64f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level, Ipp64fc value) {
			return ippsThreshold_LTVal_64fc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_LTVal_16s(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16s level, Ipp16sc value) {
			return ippsThreshold_LTVal_16sc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_LTVal_32f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level, Ipp32fc value) {
			return ippsThreshold_LTVal_32fc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_LTVal_64f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTVal(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level, Ipp64fc value) {
			return ippsThreshold_LTVal_64fc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal_16s(Ipp16s* pSrcDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_GTVal_16s_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal_16sc(Ipp16sc* pSrcDst, int len, Ipp16s level, Ipp16sc value) {
			return ippsThreshold_GTVal_16sc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal_32f(Ipp32f* pSrcDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_GTVal_32f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level, Ipp32fc value) {
			return ippsThreshold_GTVal_32fc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal_64f(Ipp64f* pSrcDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_GTVal_64f_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level, Ipp64fc value) {
			return ippsThreshold_GTVal_64fc_I(pSrcDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s level, Ipp16s value) {
			return ippsThreshold_GTVal_16s(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp16s level, Ipp16sc value) {
			return ippsThreshold_GTVal_16sc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level, Ipp32f value) {
			return ippsThreshold_GTVal_32f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level, Ipp32fc value) {
			return ippsThreshold_GTVal_32fc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level, Ipp64f value) {
			return ippsThreshold_GTVal_64f(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_GTVal(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level, Ipp64fc value) {
			return ippsThreshold_GTVal_64fc(pSrc, pDst, len, level, value);
		}

		inline auto Threshold_LTValGTVal_16s(Ipp16s* pSrcDst, int len, Ipp16s levelLT, Ipp16s valueLT, Ipp16s levelGT, Ipp16s valueGT) {
			return ippsThreshold_LTValGTVal_16s_I(pSrcDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp16s levelLT, Ipp16s valueLT, Ipp16s levelGT, Ipp16s valueGT) {
			return ippsThreshold_LTValGTVal_16s(pSrc, pDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal_32s(Ipp32s* pSrcDst, int len, Ipp32s levelLT, Ipp32s valueLT, Ipp32s levelGT, Ipp32s valueGT) {
			return ippsThreshold_LTValGTVal_32s_I(pSrcDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal(const Ipp32s* pSrc, Ipp32s* pDst, int len, Ipp32s levelLT, Ipp32s valueLT, Ipp32s levelGT, Ipp32s valueGT) {
			return ippsThreshold_LTValGTVal_32s(pSrc, pDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal_32f(Ipp32f* pSrcDst, int len, Ipp32f levelLT, Ipp32f valueLT, Ipp32f levelGT, Ipp32f valueGT) {
			return ippsThreshold_LTValGTVal_32f_I(pSrcDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f levelLT, Ipp32f valueLT, Ipp32f levelGT, Ipp32f valueGT) {
			return ippsThreshold_LTValGTVal_32f(pSrc, pDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal_64f(Ipp64f* pSrcDst, int len, Ipp64f levelLT, Ipp64f valueLT, Ipp64f levelGT, Ipp64f valueGT) {
			return ippsThreshold_LTValGTVal_64f_I(pSrcDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTValGTVal(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f levelLT, Ipp64f valueLT, Ipp64f levelGT, Ipp64f valueGT) {
			return ippsThreshold_LTValGTVal_64f(pSrc, pDst, len, levelLT, valueLT, levelGT, valueGT);
		}

		inline auto Threshold_LTInv_32f(Ipp32f* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_LTInv_32f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTInv_32fc(Ipp32fc* pSrcDst, int len, Ipp32f level) {
			return ippsThreshold_LTInv_32fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTInv_64f(Ipp64f* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_LTInv_64f_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTInv_64fc(Ipp64fc* pSrcDst, int len, Ipp64f level) {
			return ippsThreshold_LTInv_64fc_I(pSrcDst, len, level);
		}

		inline auto Threshold_LTInv(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f level) {
			return ippsThreshold_LTInv_32f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTInv(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f level) {
			return ippsThreshold_LTInv_32fc(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTInv(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f level) {
			return ippsThreshold_LTInv_64f(pSrc, pDst, len, level);
		}

		inline auto Threshold_LTInv(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f level) {
			return ippsThreshold_LTInv_64fc(pSrc, pDst, len, level);
		}

		inline auto CartToPolar(const Ipp32fc* pSrc, Ipp32f* pDstMagn, Ipp32f* pDstPhase, int len) {
			return ippsCartToPolar_32fc(pSrc, pDstMagn, pDstPhase, len);
		}

		inline auto CartToPolar(const Ipp64fc* pSrc, Ipp64f* pDstMagn, Ipp64f* pDstPhase, int len) {
			return ippsCartToPolar_64fc(pSrc, pDstMagn, pDstPhase, len);
		}

		inline auto CartToPolar(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDstMagn, Ipp32f* pDstPhase, int len) {
			return ippsCartToPolar_32f(pSrcRe, pSrcIm, pDstMagn, pDstPhase, len);
		}

		inline auto CartToPolar(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDstMagn, Ipp64f* pDstPhase, int len) {
			return ippsCartToPolar_64f(pSrcRe, pSrcIm, pDstMagn, pDstPhase, len);
		}

		inline auto CartToPolar(const Ipp16sc* pSrc, Ipp16s* pDstMagn, Ipp16s* pDstPhase, int len, int magnScaleFactor, int phaseScaleFactor) {
			return ippsCartToPolar_16sc_Sfs(pSrc, pDstMagn, pDstPhase, len, magnScaleFactor, phaseScaleFactor);
		}

		inline auto PolarToCart(const Ipp32f* pSrcMagn, const Ipp32f* pSrcPhase, Ipp32fc* pDst, int len) {
			return ippsPolarToCart_32fc(pSrcMagn, pSrcPhase, pDst, len);
		}

		inline auto PolarToCart(const Ipp64f* pSrcMagn, const Ipp64f* pSrcPhase, Ipp64fc* pDst, int len) {
			return ippsPolarToCart_64fc(pSrcMagn, pSrcPhase, pDst, len);
		}

		inline auto PolarToCart(const Ipp32f* pSrcMagn, const Ipp32f* pSrcPhase, Ipp32f* pDstRe, Ipp32f* pDstIm, int len) {
			return ippsPolarToCart_32f(pSrcMagn, pSrcPhase, pDstRe, pDstIm, len);
		}

		inline auto PolarToCart(const Ipp64f* pSrcMagn, const Ipp64f* pSrcPhase, Ipp64f* pDstRe, Ipp64f* pDstIm, int len) {
			return ippsPolarToCart_64f(pSrcMagn, pSrcPhase, pDstRe, pDstIm, len);
		}

		inline auto PolarToCart(const Ipp16s* pSrcMagn, const Ipp16s* pSrcPhase, Ipp16sc* pDst, int len, int magnScaleFactor, int phaseScaleFactor) {
			return ippsPolarToCart_16sc_Sfs(pSrcMagn, pSrcPhase, pDst, len, magnScaleFactor, phaseScaleFactor);
		}

		inline auto Flip(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
			return ippsFlip_8u(pSrc, pDst, len);
		}

		inline auto Flip_8u(Ipp8u* pSrcDst, int len) {
			return ippsFlip_8u_I(pSrcDst, len);
		}

		inline auto Flip(const Ipp16u* pSrc, Ipp16u* pDst, int len) {
			return ippsFlip_16u(pSrc, pDst, len);
		}

		inline auto Flip_16u(Ipp16u* pSrcDst, int len) {
			return ippsFlip_16u_I(pSrcDst, len);
		}

		inline auto Flip(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsFlip_32f(pSrc, pDst, len);
		}

		inline auto Flip_32f(Ipp32f* pSrcDst, int len) {
			return ippsFlip_32f_I(pSrcDst, len);
		}

		inline auto Flip(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsFlip_32fc(pSrc, pDst, len);
		}

		inline auto Flip_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsFlip_32fc_I(pSrcDst, len);
		}

		inline auto Flip(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsFlip_64f(pSrc, pDst, len);
		}

		inline auto Flip_64f(Ipp64f* pSrcDst, int len) {
			return ippsFlip_64f_I(pSrcDst, len);
		}

		inline auto Flip(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsFlip_64fc(pSrc, pDst, len);
		}

		inline auto Flip_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsFlip_64fc_I(pSrcDst, len);
		}

		inline auto FindNearestOne(Ipp16u inpVal, Ipp16u* pOutVal, int* pOutIndex, const Ipp16u* pTable, int tblLen) {
			return ippsFindNearestOne_16u(inpVal, pOutVal, pOutIndex, pTable, tblLen);
		}

		inline auto FindNearest(const Ipp16u* pVals, Ipp16u* pOutVals, int* pOutIndexes, int len, const Ipp16u* pTable, int tblLen) {
			return ippsFindNearest_16u(pVals, pOutVals, pOutIndexes, len, pTable, tblLen);
		}

		inline auto WinBartlett_16s(Ipp16s* pSrcDst, int len) {
			return ippsWinBartlett_16s_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsWinBartlett_16s(pSrc, pDst, len);
		}

		inline auto WinBartlett_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsWinBartlett_16sc_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsWinBartlett_16sc(pSrc, pDst, len);
		}

		inline auto WinBartlett_32f(Ipp32f* pSrcDst, int len) {
			return ippsWinBartlett_32f_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsWinBartlett_32f(pSrc, pDst, len);
		}

		inline auto WinBartlett_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsWinBartlett_32fc_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsWinBartlett_32fc(pSrc, pDst, len);
		}

		inline auto WinBartlett_64f(Ipp64f* pSrcDst, int len) {
			return ippsWinBartlett_64f_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsWinBartlett_64f(pSrc, pDst, len);
		}

		inline auto WinBartlett_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsWinBartlett_64fc_I(pSrcDst, len);
		}

		inline auto WinBartlett(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsWinBartlett_64fc(pSrc, pDst, len);
		}

		inline auto WinHann_16s(Ipp16s* pSrcDst, int len) {
			return ippsWinHann_16s_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsWinHann_16s(pSrc, pDst, len);
		}

		inline auto WinHann_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsWinHann_16sc_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsWinHann_16sc(pSrc, pDst, len);
		}

		inline auto WinHann_32f(Ipp32f* pSrcDst, int len) {
			return ippsWinHann_32f_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsWinHann_32f(pSrc, pDst, len);
		}

		inline auto WinHann_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsWinHann_32fc_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsWinHann_32fc(pSrc, pDst, len);
		}

		inline auto WinHann_64f(Ipp64f* pSrcDst, int len) {
			return ippsWinHann_64f_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsWinHann_64f(pSrc, pDst, len);
		}

		inline auto WinHann_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsWinHann_64fc_I(pSrcDst, len);
		}

		inline auto WinHann(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsWinHann_64fc(pSrc, pDst, len);
		}

		inline auto WinHamming_16s(Ipp16s* pSrcDst, int len) {
			return ippsWinHamming_16s_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsWinHamming_16s(pSrc, pDst, len);
		}

		inline auto WinHamming_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsWinHamming_16sc_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsWinHamming_16sc(pSrc, pDst, len);
		}

		inline auto WinHamming_32f(Ipp32f* pSrcDst, int len) {
			return ippsWinHamming_32f_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsWinHamming_32f(pSrc, pDst, len);
		}

		inline auto WinHamming_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsWinHamming_32fc_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsWinHamming_32fc(pSrc, pDst, len);
		}

		inline auto WinHamming_64f(Ipp64f* pSrcDst, int len) {
			return ippsWinHamming_64f_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsWinHamming_64f(pSrc, pDst, len);
		}

		inline auto WinHamming_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsWinHamming_64fc_I(pSrcDst, len);
		}

		inline auto WinHamming(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsWinHamming_64fc(pSrc, pDst, len);
		}

		inline auto WinBlackman_16s(Ipp16s* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_16s_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_16s(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackman_16sc(Ipp16sc* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_16sc_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_16sc(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackman_32f(Ipp32f* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_32f_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_32f(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackman_32fc(Ipp32fc* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_32fc_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f alpha) {
			return ippsWinBlackman_32fc(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackman_64f(Ipp64f* pSrcDst, int len, Ipp64f alpha) {
			return ippsWinBlackman_64f_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f alpha) {
			return ippsWinBlackman_64f(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackman_64fc(Ipp64fc* pSrcDst, int len, Ipp64f alpha) {
			return ippsWinBlackman_64fc_I(pSrcDst, len, alpha);
		}

		inline auto WinBlackman(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f alpha) {
			return ippsWinBlackman_64fc(pSrc, pDst, len, alpha);
		}

		inline auto WinBlackmanStd_16s(Ipp16s* pSrcDst, int len) {
			return ippsWinBlackmanStd_16s_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsWinBlackmanStd_16s(pSrc, pDst, len);
		}

		inline auto WinBlackmanStd_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsWinBlackmanStd_16sc_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsWinBlackmanStd_16sc(pSrc, pDst, len);
		}

		inline auto WinBlackmanStd_32f(Ipp32f* pSrcDst, int len) {
			return ippsWinBlackmanStd_32f_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsWinBlackmanStd_32f(pSrc, pDst, len);
		}

		inline auto WinBlackmanStd_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsWinBlackmanStd_32fc_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsWinBlackmanStd_32fc(pSrc, pDst, len);
		}

		inline auto WinBlackmanStd_64f(Ipp64f* pSrcDst, int len) {
			return ippsWinBlackmanStd_64f_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsWinBlackmanStd_64f(pSrc, pDst, len);
		}

		inline auto WinBlackmanStd_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsWinBlackmanStd_64fc_I(pSrcDst, len);
		}

		inline auto WinBlackmanStd(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsWinBlackmanStd_64fc(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_16s(Ipp16s* pSrcDst, int len) {
			return ippsWinBlackmanOpt_16s_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
			return ippsWinBlackmanOpt_16s(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_16sc(Ipp16sc* pSrcDst, int len) {
			return ippsWinBlackmanOpt_16sc_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
			return ippsWinBlackmanOpt_16sc(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_32f(Ipp32f* pSrcDst, int len) {
			return ippsWinBlackmanOpt_32f_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
			return ippsWinBlackmanOpt_32f(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_32fc(Ipp32fc* pSrcDst, int len) {
			return ippsWinBlackmanOpt_32fc_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
			return ippsWinBlackmanOpt_32fc(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_64f(Ipp64f* pSrcDst, int len) {
			return ippsWinBlackmanOpt_64f_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
			return ippsWinBlackmanOpt_64f(pSrc, pDst, len);
		}

		inline auto WinBlackmanOpt_64fc(Ipp64fc* pSrcDst, int len) {
			return ippsWinBlackmanOpt_64fc_I(pSrcDst, len);
		}

		inline auto WinBlackmanOpt(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
			return ippsWinBlackmanOpt_64fc(pSrc, pDst, len);
		}

		inline auto WinKaiser_16s(Ipp16s* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_16s_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp16s* pSrc, Ipp16s* pDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_16s(pSrc, pDst, len, alpha);
		}

		inline auto WinKaiser_16sc(Ipp16sc* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_16sc_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_16sc(pSrc, pDst, len, alpha);
		}

		inline auto WinKaiser_32f(Ipp32f* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_32f_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp32f* pSrc, Ipp32f* pDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_32f(pSrc, pDst, len, alpha);
		}

		inline auto WinKaiser_32fc(Ipp32fc* pSrcDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_32fc_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, Ipp32f alpha) {
			return ippsWinKaiser_32fc(pSrc, pDst, len, alpha);
		}

		inline auto WinKaiser_64f(Ipp64f* pSrcDst, int len, Ipp64f alpha) {
			return ippsWinKaiser_64f_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp64f* pSrc, Ipp64f* pDst, int len, Ipp64f alpha) {
			return ippsWinKaiser_64f(pSrc, pDst, len, alpha);
		}

		inline auto WinKaiser_64fc(Ipp64fc* pSrcDst, int len, Ipp64f alpha) {
			return ippsWinKaiser_64fc_I(pSrcDst, len, alpha);
		}

		inline auto WinKaiser(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, Ipp64f alpha) {
			return ippsWinKaiser_64fc(pSrc, pDst, len, alpha);
		}

		inline auto Sum(const Ipp16s* pSrc, int len, Ipp16s* pSum, int scaleFactor) {
			return ippsSum_16s_Sfs(pSrc, len, pSum, scaleFactor);
		}

		inline auto Sum(const Ipp16sc* pSrc, int len, Ipp16sc* pSum, int scaleFactor) {
			return ippsSum_16sc_Sfs(pSrc, len, pSum, scaleFactor);
		}

		inline auto Sum(const Ipp16s* pSrc, int len, Ipp32s* pSum, int scaleFactor) {
			return ippsSum_16s32s_Sfs(pSrc, len, pSum, scaleFactor);
		}

		inline auto Sum(const Ipp16sc* pSrc, int len, Ipp32sc* pSum, int scaleFactor) {
			return ippsSum_16sc32sc_Sfs(pSrc, len, pSum, scaleFactor);
		}

		inline auto Sum(const Ipp32s* pSrc, int len, Ipp32s* pSum, int scaleFactor) {
			return ippsSum_32s_Sfs(pSrc, len, pSum, scaleFactor);
		}

		inline auto Sum(const Ipp32f* pSrc, int len, Ipp32f* pSum, IppHintAlgorithm hint) {
			return ippsSum_32f(pSrc, len, pSum, hint);
		}

		inline auto Sum(const Ipp32fc* pSrc, int len, Ipp32fc* pSum, IppHintAlgorithm hint) {
			return ippsSum_32fc(pSrc, len, pSum, hint);
		}

		inline auto Sum(const Ipp64f* pSrc, int len, Ipp64f* pSum) {
			return ippsSum_64f(pSrc, len, pSum);
		}

		inline auto Sum(const Ipp64fc* pSrc, int len, Ipp64fc* pSum) {
			return ippsSum_64fc(pSrc, len, pSum);
		}

		inline auto Min(const Ipp16s* pSrc, int len, Ipp16s* pMin) {
			return ippsMin_16s(pSrc, len, pMin);
		}

		inline auto Min(const Ipp32s* pSrc, int len, Ipp32s* pMin) {
			return ippsMin_32s(pSrc, len, pMin);
		}

		inline auto Min(const Ipp32f* pSrc, int len, Ipp32f* pMin) {
			return ippsMin_32f(pSrc, len, pMin);
		}

		inline auto Min(const Ipp64f* pSrc, int len, Ipp64f* pMin) {
			return ippsMin_64f(pSrc, len, pMin);
		}

		inline auto Max(const Ipp16s* pSrc, int len, Ipp16s* pMax) {
			return ippsMax_16s(pSrc, len, pMax);
		}

		inline auto Max(const Ipp32s* pSrc, int len, Ipp32s* pMax) {
			return ippsMax_32s(pSrc, len, pMax);
		}

		inline auto Max(const Ipp32f* pSrc, int len, Ipp32f* pMax) {
			return ippsMax_32f(pSrc, len, pMax);
		}

		inline auto Max(const Ipp64f* pSrc, int len, Ipp64f* pMax) {
			return ippsMax_64f(pSrc, len, pMax);
		}

		inline auto MinMax(const Ipp8u* pSrc, int len, Ipp8u* pMin, Ipp8u* pMax) {
			return ippsMinMax_8u(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp16u* pSrc, int len, Ipp16u* pMin, Ipp16u* pMax) {
			return ippsMinMax_16u(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp16s* pSrc, int len, Ipp16s* pMin, Ipp16s* pMax) {
			return ippsMinMax_16s(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp32u* pSrc, int len, Ipp32u* pMin, Ipp32u* pMax) {
			return ippsMinMax_32u(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp32s* pSrc, int len, Ipp32s* pMin, Ipp32s* pMax) {
			return ippsMinMax_32s(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp32f* pSrc, int len, Ipp32f* pMin, Ipp32f* pMax) {
			return ippsMinMax_32f(pSrc, len, pMin, pMax);
		}

		inline auto MinMax(const Ipp64f* pSrc, int len, Ipp64f* pMin, Ipp64f* pMax) {
			return ippsMinMax_64f(pSrc, len, pMin, pMax);
		}

		inline auto MinAbs(const Ipp16s* pSrc, int len, Ipp16s* pMinAbs) {
			return ippsMinAbs_16s(pSrc, len, pMinAbs);
		}

		inline auto MinAbs(const Ipp32s* pSrc, int len, Ipp32s* pMinAbs) {
			return ippsMinAbs_32s(pSrc, len, pMinAbs);
		}

		inline auto MinAbs(const Ipp32f* pSrc, int len, Ipp32f* pMinAbs) {
			return ippsMinAbs_32f(pSrc, len, pMinAbs);
		}

		inline auto MinAbs(const Ipp64f* pSrc, int len, Ipp64f* pMinAbs) {
			return ippsMinAbs_64f(pSrc, len, pMinAbs);
		}

		inline auto MaxAbs(const Ipp16s* pSrc, int len, Ipp16s* pMaxAbs) {
			return ippsMaxAbs_16s(pSrc, len, pMaxAbs);
		}

		inline auto MaxAbs(const Ipp32s* pSrc, int len, Ipp32s* pMaxAbs) {
			return ippsMaxAbs_32s(pSrc, len, pMaxAbs);
		}

		inline auto MaxAbs(const Ipp32f* pSrc, int len, Ipp32f* pMaxAbs) {
			return ippsMaxAbs_32f(pSrc, len, pMaxAbs);
		}

		inline auto MaxAbs(const Ipp64f* pSrc, int len, Ipp64f* pMaxAbs) {
			return ippsMaxAbs_64f(pSrc, len, pMaxAbs);
		}

		inline auto MinIndx(const Ipp16s* pSrc, int len, Ipp16s* pMin, int* pIndx) {
			return ippsMinIndx_16s(pSrc, len, pMin, pIndx);
		}

		inline auto MinIndx(const Ipp32s* pSrc, int len, Ipp32s* pMin, int* pIndx) {
			return ippsMinIndx_32s(pSrc, len, pMin, pIndx);
		}

		inline auto MinIndx(const Ipp32f* pSrc, int len, Ipp32f* pMin, int* pIndx) {
			return ippsMinIndx_32f(pSrc, len, pMin, pIndx);
		}

		inline auto MinIndx(const Ipp64f* pSrc, int len, Ipp64f* pMin, int* pIndx) {
			return ippsMinIndx_64f(pSrc, len, pMin, pIndx);
		}

		inline auto MaxIndx(const Ipp16s* pSrc, int len, Ipp16s* pMax, int* pIndx) {
			return ippsMaxIndx_16s(pSrc, len, pMax, pIndx);
		}

		inline auto MaxIndx(const Ipp32s* pSrc, int len, Ipp32s* pMax, int* pIndx) {
			return ippsMaxIndx_32s(pSrc, len, pMax, pIndx);
		}

		inline auto MaxIndx(const Ipp32f* pSrc, int len, Ipp32f* pMax, int* pIndx) {
			return ippsMaxIndx_32f(pSrc, len, pMax, pIndx);
		}

		inline auto MaxIndx(const Ipp64f* pSrc, int len, Ipp64f* pMax, int* pIndx) {
			return ippsMaxIndx_64f(pSrc, len, pMax, pIndx);
		}

		inline auto MinMaxIndx(const Ipp8u* pSrc, int len, Ipp8u* pMin, int* pMinIndx, Ipp8u* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_8u(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp16u* pSrc, int len, Ipp16u* pMin, int* pMinIndx, Ipp16u* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_16u(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp16s* pSrc, int len, Ipp16s* pMin, int* pMinIndx, Ipp16s* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_16s(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp32u* pSrc, int len, Ipp32u* pMin, int* pMinIndx, Ipp32u* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_32u(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp32s* pSrc, int len, Ipp32s* pMin, int* pMinIndx, Ipp32s* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_32s(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp32f* pSrc, int len, Ipp32f* pMin, int* pMinIndx, Ipp32f* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_32f(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinMaxIndx(const Ipp64f* pSrc, int len, Ipp64f* pMin, int* pMinIndx, Ipp64f* pMax, int* pMaxIndx) {
			return ippsMinMaxIndx_64f(pSrc, len, pMin, pMinIndx, pMax, pMaxIndx);
		}

		inline auto MinAbsIndx(const Ipp16s* pSrc, int len, Ipp16s* pMinAbs, int* pIndx) {
			return ippsMinAbsIndx_16s(pSrc, len, pMinAbs, pIndx);
		}

		inline auto MinAbsIndx(const Ipp32s* pSrc, int len, Ipp32s* pMinAbs, int* pIndx) {
			return ippsMinAbsIndx_32s(pSrc, len, pMinAbs, pIndx);
		}

		inline auto MinAbsIndx(const Ipp32f* pSrc, int len, Ipp32f* pMinAbs, int* pIndx) {
			return ippsMinAbsIndx_32f(pSrc, len, pMinAbs, pIndx);
		}

		inline auto MinAbsIndx(const Ipp64f* pSrc, int len, Ipp64f* pMinAbs, int* pIndx) {
			return ippsMinAbsIndx_64f(pSrc, len, pMinAbs, pIndx);
		}

		inline auto MaxAbsIndx(const Ipp16s* pSrc, int len, Ipp16s* pMaxAbs, int* pIndx) {
			return ippsMaxAbsIndx_16s(pSrc, len, pMaxAbs, pIndx);
		}

		inline auto MaxAbsIndx(const Ipp32s* pSrc, int len, Ipp32s* pMaxAbs, int* pIndx) {
			return ippsMaxAbsIndx_32s(pSrc, len, pMaxAbs, pIndx);
		}

		inline auto MaxAbsIndx(const Ipp32f* pSrc, int len, Ipp32f* pMaxAbs, int* pIndx) {
			return ippsMaxAbsIndx_32f(pSrc, len, pMaxAbs, pIndx);
		}

		inline auto MaxAbsIndx(const Ipp64f* pSrc, int len, Ipp64f* pMaxAbs, int* pIndx) {
			return ippsMaxAbsIndx_64f(pSrc, len, pMaxAbs, pIndx);
		}

		inline auto Mean(const Ipp16s* pSrc, int len, Ipp16s* pMean, int scaleFactor) {
			return ippsMean_16s_Sfs(pSrc, len, pMean, scaleFactor);
		}

		inline auto Mean(const Ipp16sc* pSrc, int len, Ipp16sc* pMean, int scaleFactor) {
			return ippsMean_16sc_Sfs(pSrc, len, pMean, scaleFactor);
		}

		inline auto Mean(const Ipp32s* pSrc, int len, Ipp32s* pMean, int scaleFactor) {
			return ippsMean_32s_Sfs(pSrc, len, pMean, scaleFactor);
		}

		inline auto Mean(const Ipp32f* pSrc, int len, Ipp32f* pMean, IppHintAlgorithm hint) {
			return ippsMean_32f(pSrc, len, pMean, hint);
		}

		inline auto Mean(const Ipp32fc* pSrc, int len, Ipp32fc* pMean, IppHintAlgorithm hint) {
			return ippsMean_32fc(pSrc, len, pMean, hint);
		}

		inline auto Mean(const Ipp64f* pSrc, int len, Ipp64f* pMean) {
			return ippsMean_64f(pSrc, len, pMean);
		}

		inline auto Mean(const Ipp64fc* pSrc, int len, Ipp64fc* pMean) {
			return ippsMean_64fc(pSrc, len, pMean);
		}

		inline auto StdDev(const Ipp16s* pSrc, int len, Ipp16s* pStdDev, int scaleFactor) {
			return ippsStdDev_16s_Sfs(pSrc, len, pStdDev, scaleFactor);
		}

		inline auto StdDev(const Ipp16s* pSrc, int len, Ipp32s* pStdDev, int scaleFactor) {
			return ippsStdDev_16s32s_Sfs(pSrc, len, pStdDev, scaleFactor);
		}

		inline auto StdDev(const Ipp32f* pSrc, int len, Ipp32f* pStdDev, IppHintAlgorithm hint) {
			return ippsStdDev_32f(pSrc, len, pStdDev, hint);
		}

		inline auto StdDev(const Ipp64f* pSrc, int len, Ipp64f* pStdDev) {
			return ippsStdDev_64f(pSrc, len, pStdDev);
		}

		inline auto MeanStdDev(const Ipp16s* pSrc, int len, Ipp16s* pMean, Ipp16s* pStdDev, int scaleFactor) {
			return ippsMeanStdDev_16s_Sfs(pSrc, len, pMean, pStdDev, scaleFactor);
		}

		inline auto MeanStdDev(const Ipp16s* pSrc, int len, Ipp32s* pMean, Ipp32s* pStdDev, int scaleFactor) {
			return ippsMeanStdDev_16s32s_Sfs(pSrc, len, pMean, pStdDev, scaleFactor);
		}

		inline auto MeanStdDev(const Ipp32f* pSrc, int len, Ipp32f* pMean, Ipp32f* pStdDev, IppHintAlgorithm hint) {
			return ippsMeanStdDev_32f(pSrc, len, pMean, pStdDev, hint);
		}

		inline auto MeanStdDev(const Ipp64f* pSrc, int len, Ipp64f* pMean, Ipp64f* pStdDev) {
			return ippsMeanStdDev_64f(pSrc, len, pMean, pStdDev);
		}

		inline auto Norm_Inf(const Ipp16s* pSrc, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNorm_Inf_16s32s_Sfs(pSrc, len, pNorm, scaleFactor);
		}

		inline auto Norm_Inf(const Ipp16s* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_Inf_16s32f(pSrc, len, pNorm);
		}

		inline auto Norm_Inf(const Ipp32f* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_Inf_32f(pSrc, len, pNorm);
		}

		inline auto Norm_Inf(const Ipp32fc* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_Inf_32fc32f(pSrc, len, pNorm);
		}

		inline auto Norm_Inf(const Ipp64f* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_Inf_64f(pSrc, len, pNorm);
		}

		inline auto Norm_Inf(const Ipp64fc* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_Inf_64fc64f(pSrc, len, pNorm);
		}

		inline auto Norm_L1(const Ipp16s* pSrc, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNorm_L1_16s32s_Sfs(pSrc, len, pNorm, scaleFactor);
		}

		inline auto Norm_L1(const Ipp16s* pSrc, int len, Ipp64s* pNorm, int scaleFactor) {
			return ippsNorm_L1_16s64s_Sfs(pSrc, len, pNorm, scaleFactor);
		}

		inline auto Norm_L1(const Ipp16s* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_L1_16s32f(pSrc, len, pNorm);
		}

		inline auto Norm_L1(const Ipp32f* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_L1_32f(pSrc, len, pNorm);
		}

		inline auto Norm_L1(const Ipp32fc* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L1_32fc64f(pSrc, len, pNorm);
		}

		inline auto Norm_L1(const Ipp64f* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L1_64f(pSrc, len, pNorm);
		}

		inline auto Norm_L1(const Ipp64fc* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L1_64fc64f(pSrc, len, pNorm);
		}

		inline auto Norm_L2(const Ipp16s* pSrc, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNorm_L2_16s32s_Sfs(pSrc, len, pNorm, scaleFactor);
		}

		inline auto Norm_L2(const Ipp16s* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_L2_16s32f(pSrc, len, pNorm);
		}

		inline auto Norm_L2(const Ipp32f* pSrc, int len, Ipp32f* pNorm) {
			return ippsNorm_L2_32f(pSrc, len, pNorm);
		}

		inline auto Norm_L2(const Ipp32fc* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L2_32fc64f(pSrc, len, pNorm);
		}

		inline auto Norm_L2(const Ipp64f* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L2_64f(pSrc, len, pNorm);
		}

		inline auto Norm_L2(const Ipp64fc* pSrc, int len, Ipp64f* pNorm) {
			return ippsNorm_L2_64fc64f(pSrc, len, pNorm);
		}

		inline auto Norm_L2Sqr(const Ipp16s* pSrc, int len, Ipp64s* pNorm, int scaleFactor) {
			return ippsNorm_L2Sqr_16s64s_Sfs(pSrc, len, pNorm, scaleFactor);
		}

		inline auto NormDiff_Inf(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNormDiff_Inf_16s32s_Sfs(pSrc1, pSrc2, len, pNorm, scaleFactor);
		}

		inline auto NormDiff_Inf(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_Inf_16s32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_Inf(const Ipp32f* pSrc1, const Ipp32f* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_Inf_32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_Inf(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_Inf_32fc32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_Inf(const Ipp64f* pSrc1, const Ipp64f* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_Inf_64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_Inf(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_Inf_64fc64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L1(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNormDiff_L1_16s32s_Sfs(pSrc1, pSrc2, len, pNorm, scaleFactor);
		}

		inline auto NormDiff_L1(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp64s* pNorm, int scaleFactor) {
			return ippsNormDiff_L1_16s64s_Sfs(pSrc1, pSrc2, len, pNorm, scaleFactor);
		}

		inline auto NormDiff_L1(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_L1_16s32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L1(const Ipp32f* pSrc1, const Ipp32f* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_L1_32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L1(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L1_32fc64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L1(const Ipp64f* pSrc1, const Ipp64f* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L1_64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L1(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L1_64fc64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32s* pNorm, int scaleFactor) {
			return ippsNormDiff_L2_16s32s_Sfs(pSrc1, pSrc2, len, pNorm, scaleFactor);
		}

		inline auto NormDiff_L2(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_L2_16s32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2(const Ipp32f* pSrc1, const Ipp32f* pSrc2, int len, Ipp32f* pNorm) {
			return ippsNormDiff_L2_32f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L2_32fc64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2(const Ipp64f* pSrc1, const Ipp64f* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L2_64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, int len, Ipp64f* pNorm) {
			return ippsNormDiff_L2_64fc64f(pSrc1, pSrc2, len, pNorm);
		}

		inline auto NormDiff_L2Sqr(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp64s* pNorm, int scaleFactor) {
			return ippsNormDiff_L2Sqr_16s64s_Sfs(pSrc1, pSrc2, len, pNorm, scaleFactor);
		}

		inline auto DotProd(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32s* pDp, int scaleFactor) {
			return ippsDotProd_16s32s_Sfs(pSrc1, pSrc2, len, pDp, scaleFactor);
		}

		inline auto DotProd(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp64s* pDp) {
			return ippsDotProd_16s64s(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp16sc* pSrc1, const Ipp16sc* pSrc2, int len, Ipp64sc* pDp) {
			return ippsDotProd_16sc64sc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp16s* pSrc1, const Ipp16sc* pSrc2, int len, Ipp64sc* pDp) {
			return ippsDotProd_16s16sc64sc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp16s* pSrc1, const Ipp32s* pSrc2, int len, Ipp32s* pDp, int scaleFactor) {
			return ippsDotProd_16s32s32s_Sfs(pSrc1, pSrc2, len, pDp, scaleFactor);
		}

		inline auto DotProd(const Ipp16s* pSrc1, const Ipp16s* pSrc2, int len, Ipp32f* pDp) {
			return ippsDotProd_16s32f(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32s* pSrc1, const Ipp32s* pSrc2, int len, Ipp32s* pDp, int scaleFactor) {
			return ippsDotProd_32s_Sfs(pSrc1, pSrc2, len, pDp, scaleFactor);
		}

		inline auto DotProd(const Ipp32f* pSrc1, const Ipp32f* pSrc2, int len, Ipp32f* pDp) {
			return ippsDotProd_32f(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, int len, Ipp32fc* pDp) {
			return ippsDotProd_32fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32f* pSrc1, const Ipp32fc* pSrc2, int len, Ipp32fc* pDp) {
			return ippsDotProd_32f32fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32f* pSrc1, const Ipp32fc* pSrc2, int len, Ipp64fc* pDp) {
			return ippsDotProd_32f32fc64fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32f* pSrc1, const Ipp32f* pSrc2, int len, Ipp64f* pDp) {
			return ippsDotProd_32f64f(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, int len, Ipp64fc* pDp) {
			return ippsDotProd_32fc64fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp64f* pSrc1, const Ipp64f* pSrc2, int len, Ipp64f* pDp) {
			return ippsDotProd_64f(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp64fc* pSrc1, const Ipp64fc* pSrc2, int len, Ipp64fc* pDp) {
			return ippsDotProd_64fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto DotProd(const Ipp64f* pSrc1, const Ipp64fc* pSrc2, int len, Ipp64fc* pDp) {
			return ippsDotProd_64f64fc(pSrc1, pSrc2, len, pDp);
		}

		inline auto MinEvery_8u(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len) {
			return ippsMinEvery_8u_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery_16u(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len) {
			return ippsMinEvery_16u_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery_16s(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len) {
			return ippsMinEvery_16s_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery_32s(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len) {
			return ippsMinEvery_32s_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMinEvery_32f_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, Ipp32u len) {
			return ippsMinEvery_64f_I(pSrc, pSrcDst, len);
		}

		inline auto MinEvery(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, Ipp32u len) {
			return ippsMinEvery_8u(pSrc1, pSrc2, pDst, len);
		}

		inline auto MinEvery(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, Ipp32u len) {
			return ippsMinEvery_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto MinEvery(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, Ipp32u len) {
			return ippsMinEvery_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MinEvery(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, Ipp32u len) {
			return ippsMinEvery_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MaxEvery_8u(const Ipp8u* pSrc, Ipp8u* pSrcDst, int len) {
			return ippsMaxEvery_8u_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery_16u(const Ipp16u* pSrc, Ipp16u* pSrcDst, int len) {
			return ippsMaxEvery_16u_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery_16s(const Ipp16s* pSrc, Ipp16s* pSrcDst, int len) {
			return ippsMaxEvery_16s_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery_32s(const Ipp32s* pSrc, Ipp32s* pSrcDst, int len) {
			return ippsMaxEvery_32s_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMaxEvery_32f_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, Ipp32u len) {
			return ippsMaxEvery_64f_I(pSrc, pSrcDst, len);
		}

		inline auto MaxEvery(const Ipp8u* pSrc1, const Ipp8u* pSrc2, Ipp8u* pDst, Ipp32u len) {
			return ippsMaxEvery_8u(pSrc1, pSrc2, pDst, len);
		}

		inline auto MaxEvery(const Ipp16u* pSrc1, const Ipp16u* pSrc2, Ipp16u* pDst, Ipp32u len) {
			return ippsMaxEvery_16u(pSrc1, pSrc2, pDst, len);
		}

		inline auto MaxEvery(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, Ipp32u len) {
			return ippsMaxEvery_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MaxEvery(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, Ipp32u len) {
			return ippsMaxEvery_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MaxOrder(const Ipp16s* pSrc, int len, int* pOrder) {
			return ippsMaxOrder_16s(pSrc, len, pOrder);
		}

		inline auto MaxOrder(const Ipp32s* pSrc, int len, int* pOrder) {
			return ippsMaxOrder_32s(pSrc, len, pOrder);
		}

		inline auto MaxOrder(const Ipp32f* pSrc, int len, int* pOrder) {
			return ippsMaxOrder_32f(pSrc, len, pOrder);
		}

		inline auto MaxOrder(const Ipp64f* pSrc, int len, int* pOrder) {
			return ippsMaxOrder_64f(pSrc, len, pOrder);
		}

		inline auto CountInRange(const Ipp32s* pSrc, int len, int* pCounts, Ipp32s lowerBound, Ipp32s upperBound) {
			return ippsCountInRange_32s(pSrc, len, pCounts, lowerBound, upperBound);
		}

		inline auto ZeroCrossing(const Ipp16s* pSrc, Ipp32u len, Ipp32f* pValZCR, IppsZCType zcType) {
			return ippsZeroCrossing_16s32f(pSrc, len, pValZCR, zcType);
		}

		inline auto ZeroCrossing(const Ipp32f* pSrc, Ipp32u len, Ipp32f* pValZCR, IppsZCType zcType) {
			return ippsZeroCrossing_32f(pSrc, len, pValZCR, zcType);
		}

		inline auto SampleUp(const Ipp16s* pSrc, int srcLen, Ipp16s* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_16s(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleUp(const Ipp16sc* pSrc, int srcLen, Ipp16sc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_16sc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleUp(const Ipp32f* pSrc, int srcLen, Ipp32f* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_32f(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleUp(const Ipp32fc* pSrc, int srcLen, Ipp32fc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_32fc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleUp(const Ipp64f* pSrc, int srcLen, Ipp64f* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_64f(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleUp(const Ipp64fc* pSrc, int srcLen, Ipp64fc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleUp_64fc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp16s* pSrc, int srcLen, Ipp16s* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_16s(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp16sc* pSrc, int srcLen, Ipp16sc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_16sc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp32f* pSrc, int srcLen, Ipp32f* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_32f(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp32fc* pSrc, int srcLen, Ipp32fc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_32fc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp64f* pSrc, int srcLen, Ipp64f* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_64f(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto SampleDown(const Ipp64fc* pSrc, int srcLen, Ipp64fc* pDst, int* pDstLen, int factor, int* pPhase) {
			return ippsSampleDown_64fc(pSrc, srcLen, pDst, pDstLen, factor, pPhase);
		}

		inline auto AutoCorrNorm(const Ipp32f* pSrc, int srcLen, Ipp32f* pDst, int dstLen, IppEnum algType, Ipp8u* pBuffer) {
			return ippsAutoCorrNorm_32f(pSrc, srcLen, pDst, dstLen, algType, pBuffer);
		}

		inline auto AutoCorrNorm(const Ipp64f* pSrc, int srcLen, Ipp64f* pDst, int dstLen, IppEnum algType, Ipp8u* pBuffer) {
			return ippsAutoCorrNorm_64f(pSrc, srcLen, pDst, dstLen, algType, pBuffer);
		}

		inline auto AutoCorrNorm(const Ipp32fc* pSrc, int srcLen, Ipp32fc* pDst, int dstLen, IppEnum algType, Ipp8u* pBuffer) {
			return ippsAutoCorrNorm_32fc(pSrc, srcLen, pDst, dstLen, algType, pBuffer);
		}

		inline auto AutoCorrNorm(const Ipp64fc* pSrc, int srcLen, Ipp64fc* pDst, int dstLen, IppEnum algType, Ipp8u* pBuffer) {
			return ippsAutoCorrNorm_64fc(pSrc, srcLen, pDst, dstLen, algType, pBuffer);
		}

		inline auto CrossCorrNorm(const Ipp32f* pSrc1, int src1Len, const Ipp32f* pSrc2, int src2Len, Ipp32f* pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer) {
			return ippsCrossCorrNorm_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, algType, pBuffer);
		}

		inline auto CrossCorrNorm(const Ipp64f* pSrc1, int src1Len, const Ipp64f* pSrc2, int src2Len, Ipp64f* pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer) {
			return ippsCrossCorrNorm_64f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, algType, pBuffer);
		}

		inline auto CrossCorrNorm(const Ipp32fc* pSrc1, int src1Len, const Ipp32fc* pSrc2, int src2Len, Ipp32fc* pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer) {
			return ippsCrossCorrNorm_32fc(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, algType, pBuffer);
		}

		inline auto CrossCorrNorm(const Ipp64fc* pSrc1, int src1Len, const Ipp64fc* pSrc2, int src2Len, Ipp64fc* pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer) {
			return ippsCrossCorrNorm_64fc(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, algType, pBuffer);
		}

		inline auto Convolve(const Ipp32f* pSrc1, int src1Len, const Ipp32f* pSrc2, int src2Len, Ipp32f* pDst, IppEnum algType, Ipp8u* pBuffer) {
			return ippsConvolve_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, algType, pBuffer);
		}

		inline auto Convolve(const Ipp64f* pSrc1, int src1Len, const Ipp64f* pSrc2, int src2Len, Ipp64f* pDst, IppEnum algType, Ipp8u* pBuffer) {
			return ippsConvolve_64f(pSrc1, src1Len, pSrc2, src2Len, pDst, algType, pBuffer);
		}

		inline auto ConvBiased(const Ipp32f* pSrc1, int src1Len, const Ipp32f* pSrc2, int src2Len, Ipp32f* pDst, int dstLen, int bias) {
			return ippsConvBiased_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, bias);
		}

		inline auto SumWindow(const Ipp8u* pSrc, Ipp32f* pDst, int len, int maskSize) {
			return ippsSumWindow_8u32f(pSrc, pDst, len, maskSize);
		}

		inline auto SumWindow(const Ipp16s* pSrc, Ipp32f* pDst, int len, int maskSize) {
			return ippsSumWindow_16s32f(pSrc, pDst, len, maskSize);
		}

		template <typename T> auto FIRSRGetSize32f(int tapsLen, int* pSpecSize, int* pBufSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FIRSRGetSize32f<Ipp32fc>(int tapsLen, int* pSpecSize, int* pBufSize) {
			return ippsFIRSRGetSize32f_32fc(tapsLen, pSpecSize, pBufSize);
		}

		inline auto FIRSRInit(const Ipp32f* pTaps, int tapsLen, IppAlgType algType, IppsFIRSpec_32f* pSpec) {
			return ippsFIRSRInit_32f(pTaps, tapsLen, algType, pSpec);
		}

		inline auto FIRSRInit(const Ipp64f* pTaps, int tapsLen, IppAlgType algType, IppsFIRSpec_64f* pSpec) {
			return ippsFIRSRInit_64f(pTaps, tapsLen, algType, pSpec);
		}

		inline auto FIRSRInit(const Ipp32fc* pTaps, int tapsLen, IppAlgType algType, IppsFIRSpec_32fc* pSpec) {
			return ippsFIRSRInit_32fc(pTaps, tapsLen, algType, pSpec);
		}

		inline auto FIRSRInit(const Ipp64fc* pTaps, int tapsLen, IppAlgType algType, IppsFIRSpec_64fc* pSpec) {
			return ippsFIRSRInit_64fc(pTaps, tapsLen, algType, pSpec);
		}

		inline auto FIRSRInit32f(const Ipp32f* pTaps, int tapsLen, IppAlgType algType, IppsFIRSpec32f_32fc* pSpec) {
			return ippsFIRSRInit32f_32fc(pTaps, tapsLen, algType, pSpec);
		}

		inline auto FIRSR(const Ipp16s* pSrc, Ipp16s* pDst, int numIters, IppsFIRSpec_32f* pSpec, const Ipp16s* pDlySrc, Ipp16s* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_16s(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR(const Ipp16sc* pSrc, Ipp16sc* pDst, int numIters, IppsFIRSpec_32fc* pSpec, const Ipp16sc* pDlySrc, Ipp16sc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_16sc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR(const Ipp32f* pSrc, Ipp32f* pDst, int numIters, IppsFIRSpec_32f* pSpec, const Ipp32f* pDlySrc, Ipp32f* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_32f(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR(const Ipp64f* pSrc, Ipp64f* pDst, int numIters, IppsFIRSpec_64f* pSpec, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_64f(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR(const Ipp32fc* pSrc, Ipp32fc* pDst, int numIters, IppsFIRSpec_32fc* pSpec, const Ipp32fc* pDlySrc, Ipp32fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_32fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR(const Ipp64fc* pSrc, Ipp64fc* pDst, int numIters, IppsFIRSpec_64fc* pSpec, const Ipp64fc* pDlySrc, Ipp64fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR_64fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRSR32f(const Ipp32fc* pSrc, Ipp32fc* pDst, int numIters, IppsFIRSpec32f_32fc* pSpec, const Ipp32fc* pDlySrc, Ipp32fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRSR32f_32fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		template <typename T> auto FIRMRGetSize32f(int tapsLen, int upFactor, int downFactor, int* pSpecSize, int* pBufSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FIRMRGetSize32f<Ipp32fc>(int tapsLen, int upFactor, int downFactor, int* pSpecSize, int* pBufSize) {
			return ippsFIRMRGetSize32f_32fc(tapsLen, upFactor, downFactor, pSpecSize, pBufSize);
		}

		inline auto FIRMRInit(const Ipp32f* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, IppsFIRSpec_32f* pSpec) {
			return ippsFIRMRInit_32f(pTaps, tapsLen, upFactor, upPhase, downFactor, downPhase, pSpec);
		}

		inline auto FIRMRInit(const Ipp64f* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, IppsFIRSpec_64f* pSpec) {
			return ippsFIRMRInit_64f(pTaps, tapsLen, upFactor, upPhase, downFactor, downPhase, pSpec);
		}

		inline auto FIRMRInit(const Ipp32fc* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, IppsFIRSpec_32fc* pSpec) {
			return ippsFIRMRInit_32fc(pTaps, tapsLen, upFactor, upPhase, downFactor, downPhase, pSpec);
		}

		inline auto FIRMRInit(const Ipp64fc* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, IppsFIRSpec_64fc* pSpec) {
			return ippsFIRMRInit_64fc(pTaps, tapsLen, upFactor, upPhase, downFactor, downPhase, pSpec);
		}

		inline auto FIRMRInit32f(const Ipp32f* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, IppsFIRSpec32f_32fc* pSpec) {
			return ippsFIRMRInit32f_32fc(pTaps, tapsLen, upFactor, upPhase, downFactor, downPhase, pSpec);
		}

		inline auto FIRMR(const Ipp16s* pSrc, Ipp16s* pDst, int numIters, IppsFIRSpec_32f* pSpec, const Ipp16s* pDlySrc, Ipp16s* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_16s(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR(const Ipp16sc* pSrc, Ipp16sc* pDst, int numIters, IppsFIRSpec_32fc* pSpec, const Ipp16sc* pDlySrc, Ipp16sc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_16sc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR(const Ipp32f* pSrc, Ipp32f* pDst, int numIters, IppsFIRSpec_32f* pSpec, const Ipp32f* pDlySrc, Ipp32f* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_32f(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR(const Ipp64f* pSrc, Ipp64f* pDst, int numIters, IppsFIRSpec_64f* pSpec, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_64f(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR(const Ipp32fc* pSrc, Ipp32fc* pDst, int numIters, IppsFIRSpec_32fc* pSpec, const Ipp32fc* pDlySrc, Ipp32fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_32fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR(const Ipp64fc* pSrc, Ipp64fc* pDst, int numIters, IppsFIRSpec_64fc* pSpec, const Ipp64fc* pDlySrc, Ipp64fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR_64fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		inline auto FIRMR32f(const Ipp32fc* pSrc, Ipp32fc* pDst, int numIters, IppsFIRSpec32f_32fc* pSpec, const Ipp32fc* pDlySrc, Ipp32fc* pDlyDst, Ipp8u* pBuf) {
			return ippsFIRMR32f_32fc(pSrc, pDst, numIters, pSpec, pDlySrc, pDlyDst, pBuf);
		}

		template <typename T> auto FIRSparseGetStateSize(int nzTapsLen, int order, int* pStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FIRSparseGetStateSize<Ipp32f>(int nzTapsLen, int order, int* pStateSize) {
			return ippsFIRSparseGetStateSize_32f(nzTapsLen, order, pStateSize);
		}

		inline auto FIRSparseInit(IppsFIRSparseState_32f** ppState, const Ipp32f* pNZTaps, const Ipp32s* pNZTapPos, int nzTapsLen, const Ipp32f* pDlyLine, Ipp8u* pBuffer) {
			return ippsFIRSparseInit_32f(ppState, pNZTaps, pNZTapPos, nzTapsLen, pDlyLine, pBuffer);
		}

		template <>
		inline auto FIRSparseGetStateSize<Ipp32fc>(int nzTapsLen, int order, int* pStateSize) {
			return ippsFIRSparseGetStateSize_32fc(nzTapsLen, order, pStateSize);
		}

		inline auto FIRSparseInit(IppsFIRSparseState_32fc** ppState, const Ipp32fc* pNZTaps, const Ipp32s* pNZTapPos, int nzTapsLen, const Ipp32fc* pDlyLine, Ipp8u* pBuffer) {
			return ippsFIRSparseInit_32fc(ppState, pNZTaps, pNZTapPos, nzTapsLen, pDlyLine, pBuffer);
		}

		inline auto FIRSparse(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsFIRSparseState_32f* pState) {
			return ippsFIRSparse_32f(pSrc, pDst, len, pState);
		}

		inline auto FIRSparse(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, IppsFIRSparseState_32fc* pState) {
			return ippsFIRSparse_32fc(pSrc, pDst, len, pState);
		}

		inline auto FIRSparseSetDlyLine(IppsFIRSparseState_32f* pState, const Ipp32f* pDlyLine) {
			return ippsFIRSparseSetDlyLine_32f(pState, pDlyLine);
		}

		inline auto FIRSparseGetDlyLine(const IppsFIRSparseState_32f* pState, Ipp32f* pDlyLine) {
			return ippsFIRSparseGetDlyLine_32f(pState, pDlyLine);
		}

		inline auto FIRSparseSetDlyLine(IppsFIRSparseState_32fc* pState, const Ipp32fc* pDlyLine) {
			return ippsFIRSparseSetDlyLine_32fc(pState, pDlyLine);
		}

		inline auto FIRSparseGetDlyLine(const IppsFIRSparseState_32fc* pState, Ipp32fc* pDlyLine) {
			return ippsFIRSparseGetDlyLine_32fc(pState, pDlyLine);
		}

		inline auto FIRGenLowpass(Ipp64f rFreq, Ipp64f* pTaps, int tapsLen, IppWinType winType, IppBool doNormal, Ipp8u* pBuffer) {
			return ippsFIRGenLowpass_64f(rFreq, pTaps, tapsLen, winType, doNormal, pBuffer);
		}

		inline auto FIRGenHighpass(Ipp64f rFreq, Ipp64f* pTaps, int tapsLen, IppWinType winType, IppBool doNormal, Ipp8u* pBuffer) {
			return ippsFIRGenHighpass_64f(rFreq, pTaps, tapsLen, winType, doNormal, pBuffer);
		}

		inline auto FIRGenBandpass(Ipp64f rLowFreq, Ipp64f rHighFreq, Ipp64f* pTaps, int tapsLen, IppWinType winType, IppBool doNormal, Ipp8u* pBuffer) {
			return ippsFIRGenBandpass_64f(rLowFreq, rHighFreq, pTaps, tapsLen, winType, doNormal, pBuffer);
		}

		inline auto FIRGenBandstop(Ipp64f rLowFreq, Ipp64f rHighFreq, Ipp64f* pTaps, int tapsLen, IppWinType winType, IppBool doNormal, Ipp8u* pBuffer) {
			return ippsFIRGenBandstop_64f(rLowFreq, rHighFreq, pTaps, tapsLen, winType, doNormal, pBuffer);
		}

		template <typename T> auto FIRLMSGetStateSize32f(int tapsLen, int dlyIndex, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FIRLMSGetStateSize32f<Ipp16s>(int tapsLen, int dlyIndex, int* pBufferSize) {
			return ippsFIRLMSGetStateSize32f_16s(tapsLen, dlyIndex, pBufferSize);
		}

		template <typename T> auto FIRLMSGetStateSize(int tapsLen, int dlyIndex, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FIRLMSGetStateSize<Ipp32f>(int tapsLen, int dlyIndex, int* pBufferSize) {
			return ippsFIRLMSGetStateSize_32f(tapsLen, dlyIndex, pBufferSize);
		}

		inline auto FIRLMSInit32f(IppsFIRLMSState32f_16s** ppState, const Ipp32f* pTaps, int tapsLen, const Ipp16s* pDlyLine, int dlyIndex, Ipp8u* pBuffer) {
			return ippsFIRLMSInit32f_16s(ppState, pTaps, tapsLen, pDlyLine, dlyIndex, pBuffer);
		}

		inline auto FIRLMSInit(IppsFIRLMSState_32f** ppState, const Ipp32f* pTaps, int tapsLen, const Ipp32f* pDlyLine, int dlyIndex, Ipp8u* pBuffer) {
			return ippsFIRLMSInit_32f(ppState, pTaps, tapsLen, pDlyLine, dlyIndex, pBuffer);
		}

		inline auto FIRLMS32f(const Ipp16s* pSrc, const Ipp16s* pRef, Ipp16s* pDst, int len, float mu, IppsFIRLMSState32f_16s* pState) {
			return ippsFIRLMS32f_16s(pSrc, pRef, pDst, len, mu, pState);
		}

		inline auto FIRLMS(const Ipp32f* pSrc, const Ipp32f* pRef, Ipp32f* pDst, int len, float mu, IppsFIRLMSState_32f* pState) {
			return ippsFIRLMS_32f(pSrc, pRef, pDst, len, mu, pState);
		}

		inline auto FIRLMSGetTaps32f(const IppsFIRLMSState32f_16s* pState, Ipp32f* pOutTaps) {
			return ippsFIRLMSGetTaps32f_16s(pState, pOutTaps);
		}

		inline auto FIRLMSGetTaps(const IppsFIRLMSState_32f* pState, Ipp32f* pOutTaps) {
			return ippsFIRLMSGetTaps_32f(pState, pOutTaps);
		}

		inline auto FIRLMSGetDlyLine32f(const IppsFIRLMSState32f_16s* pState, Ipp16s* pDlyLine, int* pDlyLineIndex) {
			return ippsFIRLMSGetDlyLine32f_16s(pState, pDlyLine, pDlyLineIndex);
		}

		inline auto FIRLMSGetDlyLine(const IppsFIRLMSState_32f* pState, Ipp32f* pDlyLine, int* pDlyLineIndex) {
			return ippsFIRLMSGetDlyLine_32f(pState, pDlyLine, pDlyLineIndex);
		}

		inline auto FIRLMSSetDlyLine32f(IppsFIRLMSState32f_16s* pState, const Ipp16s* pDlyLine, int dlyLineIndex) {
			return ippsFIRLMSSetDlyLine32f_16s(pState, pDlyLine, dlyLineIndex);
		}

		inline auto FIRLMSSetDlyLine(IppsFIRLMSState_32f* pState, const Ipp32f* pDlyLine, int dlyLineIndex) {
			return ippsFIRLMSSetDlyLine_32f(pState, pDlyLine, dlyLineIndex);
		}

		inline auto IIRGetDlyLine32f(const IppsIIRState32f_16s* pState, Ipp32f* pDlyLine) {
			return ippsIIRGetDlyLine32f_16s(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine32fc(const IppsIIRState32fc_16sc* pState, Ipp32fc* pDlyLine) {
			return ippsIIRGetDlyLine32fc_16sc(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine(const IppsIIRState_32f* pState, Ipp32f* pDlyLine) {
			return ippsIIRGetDlyLine_32f(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine(const IppsIIRState_32fc* pState, Ipp32fc* pDlyLine) {
			return ippsIIRGetDlyLine_32fc(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64f(const IppsIIRState64f_16s* pState, Ipp64f* pDlyLine) {
			return ippsIIRGetDlyLine64f_16s(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64fc(const IppsIIRState64fc_16sc* pState, Ipp64fc* pDlyLine) {
			return ippsIIRGetDlyLine64fc_16sc(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64f(const IppsIIRState64f_32s* pState, Ipp64f* pDlyLine) {
			return ippsIIRGetDlyLine64f_32s(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64fc(const IppsIIRState64fc_32sc* pState, Ipp64fc* pDlyLine) {
			return ippsIIRGetDlyLine64fc_32sc(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64f_DF1(const IppsIIRState64f_32s* pState, Ipp32s* pDlyLine) {
			return ippsIIRGetDlyLine64f_DF1_32s(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64f(const IppsIIRState64f_32f* pState, Ipp64f* pDlyLine) {
			return ippsIIRGetDlyLine64f_32f(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine64fc(const IppsIIRState64fc_32fc* pState, Ipp64fc* pDlyLine) {
			return ippsIIRGetDlyLine64fc_32fc(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine(const IppsIIRState_64f* pState, Ipp64f* pDlyLine) {
			return ippsIIRGetDlyLine_64f(pState, pDlyLine);
		}

		inline auto IIRGetDlyLine(const IppsIIRState_64fc* pState, Ipp64fc* pDlyLine) {
			return ippsIIRGetDlyLine_64fc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine32f(IppsIIRState32f_16s* pState, const Ipp32f* pDlyLine) {
			return ippsIIRSetDlyLine32f_16s(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine32fc(IppsIIRState32fc_16sc* pState, const Ipp32fc* pDlyLine) {
			return ippsIIRSetDlyLine32fc_16sc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine(IppsIIRState_32f* pState, const Ipp32f* pDlyLine) {
			return ippsIIRSetDlyLine_32f(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine(IppsIIRState_32fc* pState, const Ipp32fc* pDlyLine) {
			return ippsIIRSetDlyLine_32fc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64f(IppsIIRState64f_16s* pState, const Ipp64f* pDlyLine) {
			return ippsIIRSetDlyLine64f_16s(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64fc(IppsIIRState64fc_16sc* pState, const Ipp64fc* pDlyLine) {
			return ippsIIRSetDlyLine64fc_16sc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64f(IppsIIRState64f_32s* pState, const Ipp64f* pDlyLine) {
			return ippsIIRSetDlyLine64f_32s(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64fc(IppsIIRState64fc_32sc* pState, const Ipp64fc* pDlyLine) {
			return ippsIIRSetDlyLine64fc_32sc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64f_DF1(IppsIIRState64f_32s* pState, const Ipp32s* pDlyLine) {
			return ippsIIRSetDlyLine64f_DF1_32s(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64f(IppsIIRState64f_32f* pState, const Ipp64f* pDlyLine) {
			return ippsIIRSetDlyLine64f_32f(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine64fc(IppsIIRState64fc_32fc* pState, const Ipp64fc* pDlyLine) {
			return ippsIIRSetDlyLine64fc_32fc(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine(IppsIIRState_64f* pState, const Ipp64f* pDlyLine) {
			return ippsIIRSetDlyLine_64f(pState, pDlyLine);
		}

		inline auto IIRSetDlyLine(IppsIIRState_64fc* pState, const Ipp64fc* pDlyLine) {
			return ippsIIRSetDlyLine_64fc(pState, pDlyLine);
		}

		inline auto IIR32f(Ipp16s* pSrcDst, int len, IppsIIRState32f_16s* pState, int scaleFactor) {
			return ippsIIR32f_16s_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR32f(const Ipp16s* pSrc, Ipp16s* pDst, int len, IppsIIRState32f_16s* pState, int scaleFactor) {
			return ippsIIR32f_16s_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR32fc(Ipp16sc* pSrcDst, int len, IppsIIRState32fc_16sc* pState, int scaleFactor) {
			return ippsIIR32fc_16sc_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR32fc(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, IppsIIRState32fc_16sc* pState, int scaleFactor) {
			return ippsIIR32fc_16sc_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR_32f(Ipp32f* pSrcDst, int len, IppsIIRState_32f* pState) {
			return ippsIIR_32f_I(pSrcDst, len, pState);
		}

		inline auto IIR(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsIIRState_32f* pState) {
			return ippsIIR_32f(pSrc, pDst, len, pState);
		}

		inline auto IIR_32fc(Ipp32fc* pSrcDst, int len, IppsIIRState_32fc* pState) {
			return ippsIIR_32fc_I(pSrcDst, len, pState);
		}

		inline auto IIR(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, IppsIIRState_32fc* pState) {
			return ippsIIR_32fc(pSrc, pDst, len, pState);
		}

		inline auto IIR64f(Ipp16s* pSrcDst, int len, IppsIIRState64f_16s* pState, int scaleFactor) {
			return ippsIIR64f_16s_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR64f(const Ipp16s* pSrc, Ipp16s* pDst, int len, IppsIIRState64f_16s* pState, int scaleFactor) {
			return ippsIIR64f_16s_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR64fc(Ipp16sc* pSrcDst, int len, IppsIIRState64fc_16sc* pState, int scaleFactor) {
			return ippsIIR64fc_16sc_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR64fc(const Ipp16sc* pSrc, Ipp16sc* pDst, int len, IppsIIRState64fc_16sc* pState, int scaleFactor) {
			return ippsIIR64fc_16sc_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR64f(Ipp32s* pSrcDst, int len, IppsIIRState64f_32s* pState, int scaleFactor) {
			return ippsIIR64f_32s_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR64f(const Ipp32s* pSrc, Ipp32s* pDst, int len, IppsIIRState64f_32s* pState, int scaleFactor) {
			return ippsIIR64f_32s_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR64fc(Ipp32sc* pSrcDst, int len, IppsIIRState64fc_32sc* pState, int scaleFactor) {
			return ippsIIR64fc_32sc_ISfs(pSrcDst, len, pState, scaleFactor);
		}

		inline auto IIR64fc(const Ipp32sc* pSrc, Ipp32sc* pDst, int len, IppsIIRState64fc_32sc* pState, int scaleFactor) {
			return ippsIIR64fc_32sc_Sfs(pSrc, pDst, len, pState, scaleFactor);
		}

		inline auto IIR64f_32f(Ipp32f* pSrcDst, int len, IppsIIRState64f_32f* pState) {
			return ippsIIR64f_32f_I(pSrcDst, len, pState);
		}

		inline auto IIR64f(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsIIRState64f_32f* pState) {
			return ippsIIR64f_32f(pSrc, pDst, len, pState);
		}

		inline auto IIR64fc_32fc(Ipp32fc* pSrcDst, int len, IppsIIRState64fc_32fc* pState) {
			return ippsIIR64fc_32fc_I(pSrcDst, len, pState);
		}

		inline auto IIR64fc(const Ipp32fc* pSrc, Ipp32fc* pDst, int len, IppsIIRState64fc_32fc* pState) {
			return ippsIIR64fc_32fc(pSrc, pDst, len, pState);
		}

		inline auto IIR_64f(Ipp64f* pSrcDst, int len, IppsIIRState_64f* pState) {
			return ippsIIR_64f_I(pSrcDst, len, pState);
		}

		inline auto IIR(const Ipp64f* pSrc, Ipp64f* pDst, int len, IppsIIRState_64f* pState) {
			return ippsIIR_64f(pSrc, pDst, len, pState);
		}

		inline auto IIR_64fc(Ipp64fc* pSrcDst, int len, IppsIIRState_64fc* pState) {
			return ippsIIR_64fc_I(pSrcDst, len, pState);
		}

		inline auto IIR(const Ipp64fc* pSrc, Ipp64fc* pDst, int len, IppsIIRState_64fc* pState) {
			return ippsIIR_64fc(pSrc, pDst, len, pState);
		}

		inline auto IIR(Ipp32f** ppSrcDst, int len, int nChannels, IppsIIRState_32f** ppState) {
			return ippsIIR_32f_IP(ppSrcDst, len, nChannels, ppState);
		}

		inline auto IIR(const Ipp32f** ppSrc, Ipp32f** ppDst, int len, int nChannels, IppsIIRState_32f** ppState) {
			return ippsIIR_32f_P(ppSrc, ppDst, len, nChannels, ppState);
		}

		inline auto IIR64f(Ipp32s** ppSrcDst, int len, int nChannels, IppsIIRState64f_32s** ppState, int* pScaleFactor) {
			return ippsIIR64f_32s_IPSfs(ppSrcDst, len, nChannels, ppState, pScaleFactor);
		}

		inline auto IIR64f(const Ipp32s** ppSrc, Ipp32s** ppDst, int len, int nChannels, IppsIIRState64f_32s** ppState, int* pScaleFactor) {
			return ippsIIR64f_32s_PSfs(ppSrc, ppDst, len, nChannels, ppState, pScaleFactor);
		}

		template <typename T> auto IIRGetStateSize32f(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize32f<Ipp16s>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize32f_16s(order, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize32fc(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize32fc<Ipp16sc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize32fc_16sc(order, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize<Ipp32f>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize_32f(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize<Ipp32fc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize_32fc(order, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize64f(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize64f<Ipp16s>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64f_16s(order, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize64fc(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize64fc<Ipp16sc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_16sc(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64f<Ipp32s>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64f_32s(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64fc<Ipp32sc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_32sc(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64f<Ipp32f>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64f_32f(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64fc<Ipp32fc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_32fc(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize<Ipp64f>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize_64f(order, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize<Ipp64fc>(int order, int* pBufferSize) {
			return ippsIIRGetStateSize_64fc(order, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize32f_BiQuad(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize32f_BiQuad<Ipp16s>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize32f_BiQuad_16s(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize32fc_BiQuad(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize32fc_BiQuad<Ipp16sc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize32fc_BiQuad_16sc(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize_BiQuad(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize_BiQuad<Ipp32f>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize_BiQuad_32f(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize_BiQuad_DF1(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize_BiQuad_DF1<Ipp32f>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize_BiQuad_DF1_32f(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize_BiQuad<Ipp32fc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize_BiQuad_32fc(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize64f_BiQuad(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize64f_BiQuad<Ipp16s>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64f_BiQuad_16s(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize64fc_BiQuad(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize64fc_BiQuad<Ipp16sc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_BiQuad_16sc(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64f_BiQuad<Ipp32s>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64f_BiQuad_32s(numBq, pBufferSize);
		}

		template <typename T> auto IIRGetStateSize64f_BiQuad_DF1(int numBq, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRGetStateSize64f_BiQuad_DF1<Ipp32s>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64f_BiQuad_DF1_32s(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64fc_BiQuad<Ipp32sc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_BiQuad_32sc(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64f_BiQuad<Ipp32f>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64f_BiQuad_32f(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize64fc_BiQuad<Ipp32fc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize64fc_BiQuad_32fc(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize_BiQuad<Ipp64f>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize_BiQuad_64f(numBq, pBufferSize);
		}

		template <>
		inline auto IIRGetStateSize_BiQuad<Ipp64fc>(int numBq, int* pBufferSize) {
			return ippsIIRGetStateSize_BiQuad_64fc(numBq, pBufferSize);
		}

		inline auto IIRInit32f(IppsIIRState32f_16s** ppState, const Ipp32f* pTaps, int order, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit32f_16s(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit32fc(IppsIIRState32fc_16sc** ppState, const Ipp32fc* pTaps, int order, const Ipp32fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit32fc_16sc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit(IppsIIRState_32f** ppState, const Ipp32f* pTaps, int order, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_32f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit(IppsIIRState_32fc** ppState, const Ipp32fc* pTaps, int order, const Ipp32fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_32fc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64f(IppsIIRState64f_16s** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_16s(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc(IppsIIRState64fc_16sc** ppState, const Ipp64fc* pTaps, int order, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_16sc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64f(IppsIIRState64f_32s** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_32s(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc(IppsIIRState64fc_32sc** ppState, const Ipp64fc* pTaps, int order, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_32sc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64f(IppsIIRState64f_32f** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_32f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc(IppsIIRState64fc_32fc** ppState, const Ipp64fc* pTaps, int order, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_32fc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit(IppsIIRState_64f** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_64f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit(IppsIIRState_64fc** ppState, const Ipp64fc* pTaps, int order, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_64fc(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRInit32f_BiQuad(IppsIIRState32f_16s** ppState, const Ipp32f* pTaps, int numBq, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit32f_BiQuad_16s(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit32fc_BiQuad(IppsIIRState32fc_16sc** ppState, const Ipp32fc* pTaps, int numBq, const Ipp32fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit32fc_BiQuad_16sc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit_BiQuad(IppsIIRState_32f** ppState, const Ipp32f* pTaps, int numBq, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_BiQuad_32f(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit_BiQuad_DF1(IppsIIRState_32f** ppState, const Ipp32f* pTaps, int numBq, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_BiQuad_DF1_32f(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit_BiQuad(IppsIIRState_32fc** ppState, const Ipp32fc* pTaps, int numBq, const Ipp32fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_BiQuad_32fc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64f_BiQuad(IppsIIRState64f_16s** ppState, const Ipp64f* pTaps, int numBq, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_BiQuad_16s(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc_BiQuad(IppsIIRState64fc_16sc** ppState, const Ipp64fc* pTaps, int numBq, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_BiQuad_16sc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64f_BiQuad(IppsIIRState64f_32s** ppState, const Ipp64f* pTaps, int numBq, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_BiQuad_32s(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64f_BiQuad_DF1(IppsIIRState64f_32s** ppState, const Ipp64f* pTaps, int numBq, const Ipp32s* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_BiQuad_DF1_32s(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc_BiQuad(IppsIIRState64fc_32sc** ppState, const Ipp64fc* pTaps, int numBq, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_BiQuad_32sc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64f_BiQuad(IppsIIRState64f_32f** ppState, const Ipp64f* pTaps, int numBq, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64f_BiQuad_32f(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit64fc_BiQuad(IppsIIRState64fc_32fc** ppState, const Ipp64fc* pTaps, int numBq, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit64fc_BiQuad_32fc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit_BiQuad(IppsIIRState_64f** ppState, const Ipp64f* pTaps, int numBq, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_BiQuad_64f(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		inline auto IIRInit_BiQuad(IppsIIRState_64fc** ppState, const Ipp64fc* pTaps, int numBq, const Ipp64fc* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRInit_BiQuad_64fc(ppState, pTaps, numBq, pDlyLine, pBuf);
		}

		template <typename T> auto IIRSparseGetStateSize(int nzTapsLen1, int nzTapsLen2, int order1, int order2, int* pStateSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRSparseGetStateSize<Ipp32f>(int nzTapsLen1, int nzTapsLen2, int order1, int order2, int* pStateSize) {
			return ippsIIRSparseGetStateSize_32f(nzTapsLen1, nzTapsLen2, order1, order2, pStateSize);
		}

		inline auto IIRSparseInit(IppsIIRSparseState_32f** ppState, const Ipp32f* pNZTaps, const Ipp32s* pNZTapPos, int nzTapsLen1, int nzTapsLen2, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRSparseInit_32f(ppState, pNZTaps, pNZTapPos, nzTapsLen1, nzTapsLen2, pDlyLine, pBuf);
		}

		inline auto IIRSparse(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsIIRSparseState_32f* pState) {
			return ippsIIRSparse_32f(pSrc, pDst, len, pState);
		}

		inline auto IIRGenLowpass(Ipp64f rFreq, Ipp64f ripple, int order, Ipp64f* pTaps, IppsIIRFilterType filterType, Ipp8u* pBuffer) {
			return ippsIIRGenLowpass_64f(rFreq, ripple, order, pTaps, filterType, pBuffer);
		}

		inline auto IIRGenHighpass(Ipp64f rFreq, Ipp64f ripple, int order, Ipp64f* pTaps, IppsIIRFilterType filterType, Ipp8u* pBuffer) {
			return ippsIIRGenHighpass_64f(rFreq, ripple, order, pTaps, filterType, pBuffer);
		}

		inline auto FilterMedian_8u(Ipp8u* pSrcDst, int len, int maskSize, const Ipp8u* pDlySrc, Ipp8u* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_8u_I(pSrcDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian(const Ipp8u* pSrc, Ipp8u* pDst, int len, int maskSize, const Ipp8u* pDlySrc, Ipp8u* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_8u(pSrc, pDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian_16s(Ipp16s* pSrcDst, int len, int maskSize, const Ipp16s* pDlySrc, Ipp16s* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_16s_I(pSrcDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian(const Ipp16s* pSrc, Ipp16s* pDst, int len, int maskSize, const Ipp16s* pDlySrc, Ipp16s* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_16s(pSrc, pDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian_32s(Ipp32s* pSrcDst, int len, int maskSize, const Ipp32s* pDlySrc, Ipp32s* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_32s_I(pSrcDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian(const Ipp32s* pSrc, Ipp32s* pDst, int len, int maskSize, const Ipp32s* pDlySrc, Ipp32s* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_32s(pSrc, pDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian_32f(Ipp32f* pSrcDst, int len, int maskSize, const Ipp32f* pDlySrc, Ipp32f* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_32f_I(pSrcDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian(const Ipp32f* pSrc, Ipp32f* pDst, int len, int maskSize, const Ipp32f* pDlySrc, Ipp32f* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_32f(pSrc, pDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian_64f(Ipp64f* pSrcDst, int len, int maskSize, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_64f_I(pSrcDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto FilterMedian(const Ipp64f* pSrc, Ipp64f* pDst, int len, int maskSize, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u* pBuffer) {
			return ippsFilterMedian_64f(pSrc, pDst, len, maskSize, pDlySrc, pDlyDst, pBuffer);
		}

		inline auto ResamplePolyphase(const Ipp16s* pSrc, int len, Ipp16s* pDst, Ipp64f factor, Ipp32f norm, Ipp64f* pTime, int* pOutlen, const IppsResamplingPolyphase_16s* pSpec) {
			return ippsResamplePolyphase_16s(pSrc, len, pDst, factor, norm, pTime, pOutlen, pSpec);
		}

		inline auto ResamplePolyphase(const Ipp32f* pSrc, int len, Ipp32f* pDst, Ipp64f factor, Ipp32f norm, Ipp64f* pTime, int* pOutlen, const IppsResamplingPolyphase_32f* pSpec) {
			return ippsResamplePolyphase_32f(pSrc, len, pDst, factor, norm, pTime, pOutlen, pSpec);
		}

		inline auto ResamplePolyphaseFixed(const Ipp16s* pSrc, int len, Ipp16s* pDst, Ipp32f norm, Ipp64f* pTime, int* pOutlen, const IppsResamplingPolyphaseFixed_16s* pSpec) {
			return ippsResamplePolyphaseFixed_16s(pSrc, len, pDst, norm, pTime, pOutlen, pSpec);
		}

		inline auto ResamplePolyphaseFixed(const Ipp32f* pSrc, int len, Ipp32f* pDst, Ipp32f norm, Ipp64f* pTime, int* pOutlen, const IppsResamplingPolyphaseFixed_32f* pSpec) {
			return ippsResamplePolyphaseFixed_32f(pSrc, len, pDst, norm, pTime, pOutlen, pSpec);
		}

		template <typename T> auto ResamplePolyphaseGetSize(Ipp32f window, int nStep, int* pSize, IppHintAlgorithm hint) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto ResamplePolyphaseGetSize<Ipp16s>(Ipp32f window, int nStep, int* pSize, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseGetSize_16s(window, nStep, pSize, hint);
		}

		inline auto ResamplePolyphaseGetSize(Ipp32f window, int nStep, int* pSize, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseGetSize_32f(window, nStep, pSize, hint);
		}

		template <typename T> auto ResamplePolyphaseFixedGetSize(int inRate, int outRate, int len, int* pSize, int* pLen, int* pHeight, IppHintAlgorithm hint) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto ResamplePolyphaseFixedGetSize<Ipp16s>(int inRate, int outRate, int len, int* pSize, int* pLen, int* pHeight, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseFixedGetSize_16s(inRate, outRate, len, pSize, pLen, pHeight, hint);
		}

		template <>
		inline auto ResamplePolyphaseFixedGetSize<Ipp32f>(int inRate, int outRate, int len, int* pSize, int* pLen, int* pHeight, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseFixedGetSize_32f(inRate, outRate, len, pSize, pLen, pHeight, hint);
		}

		inline auto ResamplePolyphaseInit(Ipp32f window, int nStep, Ipp32f rollf, Ipp32f alpha, IppsResamplingPolyphase_16s* pSpec, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseInit_16s(window, nStep, rollf, alpha, pSpec, hint);
		}

		inline auto ResamplePolyphaseInit(Ipp32f window, int nStep, Ipp32f rollf, Ipp32f alpha, IppsResamplingPolyphase_32f* pSpec, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseInit_32f(window, nStep, rollf, alpha, pSpec, hint);
		}

		inline auto ResamplePolyphaseFixedInit(int inRate, int outRate, int len, Ipp32f rollf, Ipp32f alpha, IppsResamplingPolyphaseFixed_16s* pSpec, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseFixedInit_16s(inRate, outRate, len, rollf, alpha, pSpec, hint);
		}

		inline auto ResamplePolyphaseFixedInit(int inRate, int outRate, int len, Ipp32f rollf, Ipp32f alpha, IppsResamplingPolyphaseFixed_32f* pSpec, IppHintAlgorithm hint) {
			return ippsResamplePolyphaseFixedInit_32f(inRate, outRate, len, rollf, alpha, pSpec, hint);
		}

		inline auto ResamplePolyphaseSetFixedFilter(const Ipp16s* pSrc, int step, int height, IppsResamplingPolyphaseFixed_16s* pSpec) {
			return ippsResamplePolyphaseSetFixedFilter_16s(pSrc, step, height, pSpec);
		}

		inline auto ResamplePolyphaseSetFixedFilter(const Ipp32f* pSrc, int step, int height, IppsResamplingPolyphaseFixed_32f* pSpec) {
			return ippsResamplePolyphaseSetFixedFilter_32f(pSrc, step, height, pSpec);
		}

		inline auto ResamplePolyphaseGetFixedFilter(Ipp16s* pDst, int step, int height, const IppsResamplingPolyphaseFixed_16s* pSpec) {
			return ippsResamplePolyphaseGetFixedFilter_16s(pDst, step, height, pSpec);
		}

		inline auto ResamplePolyphaseGetFixedFilter(Ipp32f* pDst, int step, int height, const IppsResamplingPolyphaseFixed_32f* pSpec) {
			return ippsResamplePolyphaseGetFixedFilter_32f(pDst, step, height, pSpec);
		}

		template <typename T> auto IIRIIRGetStateSize(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRIIRGetStateSize<Ipp32f>(int order, int* pBufferSize) {
			return ippsIIRIIRGetStateSize_32f(order, pBufferSize);
		}

		template <typename T> auto IIRIIRGetStateSize64f(int order, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto IIRIIRGetStateSize64f<Ipp32f>(int order, int* pBufferSize) {
			return ippsIIRIIRGetStateSize64f_32f(order, pBufferSize);
		}

		template <>
		inline auto IIRIIRGetStateSize<Ipp64f>(int order, int* pBufferSize) {
			return ippsIIRIIRGetStateSize_64f(order, pBufferSize);
		}

		inline auto IIRIIRInit64f(IppsIIRState64f_32f** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRIIRInit64f_32f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRIIRInit(IppsIIRState_32f** ppState, const Ipp32f* pTaps, int order, const Ipp32f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRIIRInit_32f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRIIRInit(IppsIIRState_64f** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u* pBuf) {
			return ippsIIRIIRInit_64f(ppState, pTaps, order, pDlyLine, pBuf);
		}

		inline auto IIRIIR_32f(Ipp32f* pSrcDst, int len, IppsIIRState_32f* pState) {
			return ippsIIRIIR_32f_I(pSrcDst, len, pState);
		}

		inline auto IIRIIR(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsIIRState_32f* pState) {
			return ippsIIRIIR_32f(pSrc, pDst, len, pState);
		}

		inline auto IIRIIR64f_32f(Ipp32f* pSrcDst, int len, IppsIIRState64f_32f* pState) {
			return ippsIIRIIR64f_32f_I(pSrcDst, len, pState);
		}

		inline auto IIRIIR64f(const Ipp32f* pSrc, Ipp32f* pDst, int len, IppsIIRState64f_32f* pState) {
			return ippsIIRIIR64f_32f(pSrc, pDst, len, pState);
		}

		inline auto IIRIIR_64f(Ipp64f* pSrcDst, int len, IppsIIRState_64f* pState) {
			return ippsIIRIIR_64f_I(pSrcDst, len, pState);
		}

		inline auto IIRIIR(const Ipp64f* pSrc, Ipp64f* pDst, int len, IppsIIRState_64f* pState) {
			return ippsIIRIIR_64f(pSrc, pDst, len, pState);
		}

		inline auto IIRIIRGetDlyLine64f(const IppsIIRState64f_32f* pState, Ipp64f* pDlyLine) {
			return ippsIIRIIRGetDlyLine64f_32f(pState, pDlyLine);
		}

		inline auto IIRIIRSetDlyLine64f(IppsIIRState64f_32f* pState, const Ipp64f* pDlyLine) {
			return ippsIIRIIRSetDlyLine64f_32f(pState, pDlyLine);
		}

		inline auto IIRIIRGetDlyLine(const IppsIIRState_32f* pState, Ipp32f* pDlyLine) {
			return ippsIIRIIRGetDlyLine_32f(pState, pDlyLine);
		}

		inline auto IIRIIRSetDlyLine(IppsIIRState_32f* pState, const Ipp32f* pDlyLine) {
			return ippsIIRIIRSetDlyLine_32f(pState, pDlyLine);
		}

		inline auto IIRIIRGetDlyLine(const IppsIIRState_64f* pState, Ipp64f* pDlyLine) {
			return ippsIIRIIRGetDlyLine_64f(pState, pDlyLine);
		}

		inline auto IIRIIRSetDlyLine(IppsIIRState_64f* pState, const Ipp64f* pDlyLine) {
			return ippsIIRIIRSetDlyLine_64f(pState, pDlyLine);
		}

		template <typename T> auto FFTGetSize_C(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FFTGetSize_C<Ipp16fc>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_C_16fc(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto FFTGetSize_C<Ipp32f>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_C_32f(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <typename T> auto FFTGetSize_R(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto FFTGetSize_R<Ipp32f>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_R_32f(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto FFTGetSize_C<Ipp32fc>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_C_32fc(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto FFTGetSize_C<Ipp64f>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_C_64f(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto FFTGetSize_R<Ipp64f>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_R_64f(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto FFTGetSize_C<Ipp64fc>(int order, int flag, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsFFTGetSize_C_64fc(order, flag, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		inline auto FFTInit_C(IppsFFTSpec_C_16fc** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_C_16fc(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_C(IppsFFTSpec_C_32f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_C_32f(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_R(IppsFFTSpec_R_32f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_R_32f(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_C(IppsFFTSpec_C_32fc** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_C_32fc(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_C(IppsFFTSpec_C_64f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_C_64f(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_R(IppsFFTSpec_R_64f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_R_64f(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTInit_C(IppsFFTSpec_C_64fc** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsFFTInit_C_64fc(ppFFTSpec, order, flag, hint, pSpec, pSpecBuffer);
		}

		inline auto FFTFwd_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, const IppsFFTSpec_C_16fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_16fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, const IppsFFTSpec_C_16fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_16fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC(const Ipp32fc* pSrc, Ipp32fc* pDst, const IppsFFTSpec_C_32fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_32fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC(const Ipp32fc* pSrc, Ipp32fc* pDst, const IppsFFTSpec_C_32fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_32fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC(const Ipp64fc* pSrc, Ipp64fc* pDst, const IppsFFTSpec_C_64fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_64fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC(const Ipp64fc* pSrc, Ipp64fc* pDst, const IppsFFTSpec_C_64fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_64fc(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC_32fc(Ipp32fc* pSrcDst, const IppsFFTSpec_C_32fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_32fc_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC_32fc(Ipp32fc* pSrcDst, const IppsFFTSpec_C_32fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_32fc_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC_64fc(Ipp64fc* pSrcDst, const IppsFFTSpec_C_64fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_64fc_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC_64fc(Ipp64fc* pSrcDst, const IppsFFTSpec_C_64fc* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_64fc_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC_32f(Ipp32f* pSrcDstRe, Ipp32f* pSrcDstIm, const IppsFFTSpec_C_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_32f_I(pSrcDstRe, pSrcDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC_32f(Ipp32f* pSrcDstRe, Ipp32f* pSrcDstIm, const IppsFFTSpec_C_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_32f_I(pSrcDstRe, pSrcDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC_64f(Ipp64f* pSrcDstRe, Ipp64f* pSrcDstIm, const IppsFFTSpec_C_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_64f_I(pSrcDstRe, pSrcDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC_64f(Ipp64f* pSrcDstRe, Ipp64f* pSrcDstIm, const IppsFFTSpec_C_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_64f_I(pSrcDstRe, pSrcDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, const IppsFFTSpec_C_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_32f(pSrcRe, pSrcIm, pDstRe, pDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, const IppsFFTSpec_C_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_32f(pSrcRe, pSrcIm, pDstRe, pDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_CToC(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, const IppsFFTSpec_C_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_CToC_64f(pSrcRe, pSrcIm, pDstRe, pDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CToC(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, const IppsFFTSpec_C_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CToC_64f(pSrcRe, pSrcIm, pDstRe, pDstIm, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPerm_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPerm_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPack_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPack_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToCCS_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToCCS_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PermToR_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PermToR_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PackToR_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PackToR_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CCSToR_32f(Ipp32f* pSrcDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CCSToR_32f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPerm_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPerm_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPack_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPack_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToCCS_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToCCS_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PermToR_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PermToR_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PackToR_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PackToR_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CCSToR_64f(Ipp64f* pSrcDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CCSToR_64f_I(pSrcDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPerm(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPerm_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPack(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPack_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToCCS(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToCCS_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PermToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PermToR_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PackToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PackToR_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CCSToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CCSToR_32f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPerm(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPerm_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToPack(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToPack_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTFwd_RToCCS(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTFwd_RToCCS_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PermToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PermToR_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_PackToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_PackToR_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		inline auto FFTInv_CCSToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsFFTSpec_R_64f* pFFTSpec, Ipp8u* pBuffer) {
			return ippsFFTInv_CCSToR_64f(pSrc, pDst, pFFTSpec, pBuffer);
		}

		template <typename T> auto DFTGetSize_C(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DFTGetSize_C<Ipp16fc>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_C_16fc(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <>
		inline auto DFTGetSize_C<Ipp32f>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_C_32f(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <typename T> auto DFTGetSize_R(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DFTGetSize_R<Ipp32f>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_R_32f(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <>
		inline auto DFTGetSize_C<Ipp32fc>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_C_32fc(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <>
		inline auto DFTGetSize_C<Ipp64f>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_C_64f(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <>
		inline auto DFTGetSize_R<Ipp64f>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_R_64f(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		template <>
		inline auto DFTGetSize_C<Ipp64fc>(int length, int flag, IppHintAlgorithm hint, int* pSizeSpec, int* pSizeInit, int* pSizeBuf) {
			return ippsDFTGetSize_C_64fc(length, flag, hint, pSizeSpec, pSizeInit, pSizeBuf);
		}

		inline auto DFTInit_C(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_C_16fc* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_C_16fc(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_C(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_C_32f* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_C_32f(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_R(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_R_32f(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_C(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_C_32fc* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_C_32fc(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_C(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_C_64f* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_C_64f(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_R(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_R_64f(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTInit_C(int length, int flag, IppHintAlgorithm hint, IppsDFTSpec_C_64fc* pDFTSpec, Ipp8u* pMemInit) {
			return ippsDFTInit_C_64fc(length, flag, hint, pDFTSpec, pMemInit);
		}

		inline auto DFTFwd_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, const IppsDFTSpec_C_16fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_CToC_16fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, const IppsDFTSpec_C_16fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CToC_16fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_CToC(const Ipp32fc* pSrc, Ipp32fc* pDst, const IppsDFTSpec_C_32fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_CToC_32fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CToC(const Ipp32fc* pSrc, Ipp32fc* pDst, const IppsDFTSpec_C_32fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CToC_32fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_CToC(const Ipp64fc* pSrc, Ipp64fc* pDst, const IppsDFTSpec_C_64fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_CToC_64fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CToC(const Ipp64fc* pSrc, Ipp64fc* pDst, const IppsDFTSpec_C_64fc* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CToC_64fc(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_CToC(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, const IppsDFTSpec_C_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_CToC_32f(pSrcRe, pSrcIm, pDstRe, pDstIm, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CToC(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, const IppsDFTSpec_C_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CToC_32f(pSrcRe, pSrcIm, pDstRe, pDstIm, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_CToC(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, const IppsDFTSpec_C_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_CToC_64f(pSrcRe, pSrcIm, pDstRe, pDstIm, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CToC(const Ipp64f* pSrcRe, const Ipp64f* pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, const IppsDFTSpec_C_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CToC_64f(pSrcRe, pSrcIm, pDstRe, pDstIm, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_Direct_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, int length) {
			return ippsDFTFwd_Direct_CToC_16fc(pSrc, pDst, length);
		}

		inline auto DFTInv_Direct_CToC(const Ipp16fc* pSrc, Ipp16fc* pDst, int length) {
			return ippsDFTInv_Direct_CToC_16fc(pSrc, pDst, length);
		}

		inline auto DFTFwd_RToPerm(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToPerm_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_RToPack(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToPack_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_RToCCS(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToCCS_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_PermToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_PermToR_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_PackToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_PackToR_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CCSToR(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDFTSpec_R_32f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CCSToR_32f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_RToPerm(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToPerm_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_RToPack(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToPack_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTFwd_RToCCS(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTFwd_RToCCS_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_PermToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_PermToR_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_PackToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_PackToR_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto DFTInv_CCSToR(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDFTSpec_R_64f* pDFTSpec, Ipp8u* pBuffer) {
			return ippsDFTInv_CCSToR_64f(pSrc, pDst, pDFTSpec, pBuffer);
		}

		inline auto MulPack_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMulPack_32f_I(pSrc, pSrcDst, len);
		}

		inline auto MulPerm_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMulPerm_32f_I(pSrc, pSrcDst, len);
		}

		inline auto MulPack(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsMulPack_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MulPerm(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f* pDst, int len) {
			return ippsMulPerm_32f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MulPack_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsMulPack_64f_I(pSrc, pSrcDst, len);
		}

		inline auto MulPerm_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsMulPerm_64f_I(pSrc, pSrcDst, len);
		}

		inline auto MulPack(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsMulPack_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MulPerm(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f* pDst, int len) {
			return ippsMulPerm_64f(pSrc1, pSrc2, pDst, len);
		}

		inline auto MulPackConj_32f(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
			return ippsMulPackConj_32f_I(pSrc, pSrcDst, len);
		}

		inline auto MulPackConj_64f(const Ipp64f* pSrc, Ipp64f* pSrcDst, int len) {
			return ippsMulPackConj_64f_I(pSrc, pSrcDst, len);
		}

		inline auto Goertz(const Ipp16s* pSrc, int len, Ipp16sc* pVal, Ipp32f rFreq, int scaleFactor) {
			return ippsGoertz_16s_Sfs(pSrc, len, pVal, rFreq, scaleFactor);
		}

		inline auto Goertz(const Ipp16sc* pSrc, int len, Ipp16sc* pVal, Ipp32f rFreq, int scaleFactor) {
			return ippsGoertz_16sc_Sfs(pSrc, len, pVal, rFreq, scaleFactor);
		}

		inline auto Goertz(const Ipp32f* pSrc, int len, Ipp32fc* pVal, Ipp32f rFreq) {
			return ippsGoertz_32f(pSrc, len, pVal, rFreq);
		}

		inline auto Goertz(const Ipp32fc* pSrc, int len, Ipp32fc* pVal, Ipp32f rFreq) {
			return ippsGoertz_32fc(pSrc, len, pVal, rFreq);
		}

		inline auto Goertz(const Ipp64f* pSrc, int len, Ipp64fc* pVal, Ipp64f rFreq) {
			return ippsGoertz_64f(pSrc, len, pVal, rFreq);
		}

		inline auto Goertz(const Ipp64fc* pSrc, int len, Ipp64fc* pVal, Ipp64f rFreq) {
			return ippsGoertz_64fc(pSrc, len, pVal, rFreq);
		}

		template <typename T> auto DCTFwdGetSize(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DCTFwdGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsDCTFwdGetSize_32f(len, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <typename T> auto DCTInvGetSize(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto DCTInvGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsDCTInvGetSize_32f(len, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto DCTFwdGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsDCTFwdGetSize_64f(len, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		template <>
		inline auto DCTInvGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int* pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
			return ippsDCTInvGetSize_64f(len, hint, pSpecSize, pSpecBufferSize, pBufferSize);
		}

		inline auto DCTFwdInit(IppsDCTFwdSpec_32f** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsDCTFwdInit_32f(ppDCTSpec, len, hint, pSpec, pSpecBuffer);
		}

		inline auto DCTInvInit(IppsDCTInvSpec_32f** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsDCTInvInit_32f(ppDCTSpec, len, hint, pSpec, pSpecBuffer);
		}

		inline auto DCTFwdInit(IppsDCTFwdSpec_64f** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsDCTFwdInit_64f(ppDCTSpec, len, hint, pSpec, pSpecBuffer);
		}

		inline auto DCTInvInit(IppsDCTInvSpec_64f** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) {
			return ippsDCTInvInit_64f(ppDCTSpec, len, hint, pSpec, pSpecBuffer);
		}

		inline auto DCTFwd_32f(Ipp32f* pSrcDst, const IppsDCTFwdSpec_32f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTFwd_32f_I(pSrcDst, pDCTSpec, pBuffer);
		}

		inline auto DCTInv_32f(Ipp32f* pSrcDst, const IppsDCTInvSpec_32f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTInv_32f_I(pSrcDst, pDCTSpec, pBuffer);
		}

		inline auto DCTFwd(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDCTFwdSpec_32f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTFwd_32f(pSrc, pDst, pDCTSpec, pBuffer);
		}

		inline auto DCTInv(const Ipp32f* pSrc, Ipp32f* pDst, const IppsDCTInvSpec_32f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTInv_32f(pSrc, pDst, pDCTSpec, pBuffer);
		}

		inline auto DCTFwd_64f(Ipp64f* pSrcDst, const IppsDCTFwdSpec_64f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTFwd_64f_I(pSrcDst, pDCTSpec, pBuffer);
		}

		inline auto DCTInv_64f(Ipp64f* pSrcDst, const IppsDCTInvSpec_64f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTInv_64f_I(pSrcDst, pDCTSpec, pBuffer);
		}

		inline auto DCTFwd(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDCTFwdSpec_64f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTFwd_64f(pSrc, pDst, pDCTSpec, pBuffer);
		}

		inline auto DCTInv(const Ipp64f* pSrc, Ipp64f* pDst, const IppsDCTInvSpec_64f* pDCTSpec, Ipp8u* pBuffer) {
			return ippsDCTInv_64f(pSrc, pDst, pDCTSpec, pBuffer);
		}

		template <typename T> auto HilbertGetSize(int length, IppHintAlgorithm hint, int* pSpecSize, int* pBufferSize) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto HilbertGetSize<Ipp32f>(int length, IppHintAlgorithm hint, int* pSpecSize, int* pBufferSize) {
			return ippsHilbertGetSize_32f32fc(length, hint, pSpecSize, pBufferSize);
		}

		template <>
		inline auto HilbertGetSize<Ipp64f>(int length, IppHintAlgorithm hint, int* pSpecSize, int* pBufferSize) {
			return ippsHilbertGetSize_64f64fc(length, hint, pSpecSize, pBufferSize);
		}

		template <typename T> auto HilbertInit(int length, IppHintAlgorithm hint, IppsHilbertSpec* pSpec, Ipp8u* pBuffer) { static_assert(sizeof(T) == 0, "Unexpected overload"); }

		template <>
		inline auto HilbertInit<Ipp32f>(int length, IppHintAlgorithm hint, IppsHilbertSpec* pSpec, Ipp8u* pBuffer) {
			return ippsHilbertInit_32f32fc(length, hint, pSpec, pBuffer);
		}

		template <>
		inline auto HilbertInit<Ipp64f>(int length, IppHintAlgorithm hint, IppsHilbertSpec* pSpec, Ipp8u* pBuffer) {
			return ippsHilbertInit_64f64fc(length, hint, pSpec, pBuffer);
		}

		inline auto Hilbert(const Ipp32f* pSrc, Ipp32fc* pDst, IppsHilbertSpec* pSpec, Ipp8u* pBuffer) {
			return ippsHilbert_32f32fc(pSrc, pDst, pSpec, pBuffer);
		}

		inline auto Hilbert(const Ipp64f* pSrc, Ipp64fc* pDst, IppsHilbertSpec* pSpec, Ipp8u* pBuffer) {
			return ippsHilbert_64f64fc(pSrc, pDst, pSpec, pBuffer);
		}

		inline auto WTHaarFwd(const Ipp16s* pSrc, int len, Ipp16s* pDstLow, Ipp16s* pDstHigh, int scaleFactor) {
			return ippsWTHaarFwd_16s_Sfs(pSrc, len, pDstLow, pDstHigh, scaleFactor);
		}

		inline auto WTHaarFwd(const Ipp32f* pSrc, int len, Ipp32f* pDstLow, Ipp32f* pDstHigh) {
			return ippsWTHaarFwd_32f(pSrc, len, pDstLow, pDstHigh);
		}

		inline auto WTHaarFwd(const Ipp64f* pSrc, int len, Ipp64f* pDstLow, Ipp64f* pDstHigh) {
			return ippsWTHaarFwd_64f(pSrc, len, pDstLow, pDstHigh);
		}

		inline auto WTHaarInv(const Ipp16s* pSrcLow, const Ipp16s* pSrcHigh, Ipp16s* pDst, int len, int scaleFactor) {
			return ippsWTHaarInv_16s_Sfs(pSrcLow, pSrcHigh, pDst, len, scaleFactor);
		}

		inline auto WTHaarInv(const Ipp32f* pSrcLow, const Ipp32f* pSrcHigh, Ipp32f* pDst, int len) {
			return ippsWTHaarInv_32f(pSrcLow, pSrcHigh, pDst, len);
		}

		inline auto WTHaarInv(const Ipp64f* pSrcLow, const Ipp64f* pSrcHigh, Ipp64f* pDst, int len) {
			return ippsWTHaarInv_64f(pSrcLow, pSrcHigh, pDst, len);
		}

		inline auto WTFwdInit(IppsWTFwdState_8u32f* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTFwdInit_8u32f(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTFwdInit(IppsWTFwdState_16s32f* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTFwdInit_16s32f(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTFwdInit(IppsWTFwdState_16u32f* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTFwdInit_16u32f(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTFwdInit(IppsWTFwdState_32f* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTFwdInit_32f(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTFwdSetDlyLine(IppsWTFwdState_8u32f* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTFwdSetDlyLine_8u32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdSetDlyLine(IppsWTFwdState_16s32f* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTFwdSetDlyLine_16s32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdSetDlyLine(IppsWTFwdState_16u32f* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTFwdSetDlyLine_16u32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdSetDlyLine(IppsWTFwdState_32f* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTFwdSetDlyLine_32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdGetDlyLine(IppsWTFwdState_8u32f* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTFwdGetDlyLine_8u32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdGetDlyLine(IppsWTFwdState_16s32f* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTFwdGetDlyLine_16s32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdGetDlyLine(IppsWTFwdState_16u32f* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTFwdGetDlyLine_16u32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwdGetDlyLine(IppsWTFwdState_32f* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTFwdGetDlyLine_32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTFwd(const Ipp8u* pSrc, Ipp32f* pDstLow, Ipp32f* pDstHigh, int dstLen, IppsWTFwdState_8u32f* pState) {
			return ippsWTFwd_8u32f(pSrc, pDstLow, pDstHigh, dstLen, pState);
		}

		inline auto WTFwd(const Ipp16s* pSrc, Ipp32f* pDstLow, Ipp32f* pDstHigh, int dstLen, IppsWTFwdState_16s32f* pState) {
			return ippsWTFwd_16s32f(pSrc, pDstLow, pDstHigh, dstLen, pState);
		}

		inline auto WTFwd(const Ipp16u* pSrc, Ipp32f* pDstLow, Ipp32f* pDstHigh, int dstLen, IppsWTFwdState_16u32f* pState) {
			return ippsWTFwd_16u32f(pSrc, pDstLow, pDstHigh, dstLen, pState);
		}

		inline auto WTFwd(const Ipp32f* pSrc, Ipp32f* pDstLow, Ipp32f* pDstHigh, int dstLen, IppsWTFwdState_32f* pState) {
			return ippsWTFwd_32f(pSrc, pDstLow, pDstHigh, dstLen, pState);
		}

		inline auto WTInvInit(IppsWTInvState_32f8u* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTInvInit_32f8u(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTInvInit(IppsWTInvState_32f16u* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTInvInit_32f16u(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTInvInit(IppsWTInvState_32f16s* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTInvInit_32f16s(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTInvInit(IppsWTInvState_32f* pState, const Ipp32f* pTapsLow, int lenLow, int offsLow, const Ipp32f* pTapsHigh, int lenHigh, int offsHigh) {
			return ippsWTInvInit_32f(pState, pTapsLow, lenLow, offsLow, pTapsHigh, lenHigh, offsHigh);
		}

		inline auto WTInvSetDlyLine(IppsWTInvState_32f8u* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTInvSetDlyLine_32f8u(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvSetDlyLine(IppsWTInvState_32f16s* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTInvSetDlyLine_32f16s(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvSetDlyLine(IppsWTInvState_32f16u* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTInvSetDlyLine_32f16u(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvSetDlyLine(IppsWTInvState_32f* pState, const Ipp32f* pDlyLow, const Ipp32f* pDlyHigh) {
			return ippsWTInvSetDlyLine_32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvGetDlyLine(IppsWTInvState_32f8u* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTInvGetDlyLine_32f8u(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvGetDlyLine(IppsWTInvState_32f16s* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTInvGetDlyLine_32f16s(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvGetDlyLine(IppsWTInvState_32f16u* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTInvGetDlyLine_32f16u(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInvGetDlyLine(IppsWTInvState_32f* pState, Ipp32f* pDlyLow, Ipp32f* pDlyHigh) {
			return ippsWTInvGetDlyLine_32f(pState, pDlyLow, pDlyHigh);
		}

		inline auto WTInv(const Ipp32f* pSrcLow, const Ipp32f* pSrcHigh, int srcLen, Ipp8u* pDst, IppsWTInvState_32f8u* pState) {
			return ippsWTInv_32f8u(pSrcLow, pSrcHigh, srcLen, pDst, pState);
		}

		inline auto WTInv(const Ipp32f* pSrcLow, const Ipp32f* pSrcHigh, int srcLen, Ipp16s* pDst, IppsWTInvState_32f16s* pState) {
			return ippsWTInv_32f16s(pSrcLow, pSrcHigh, srcLen, pDst, pState);
		}

		inline auto WTInv(const Ipp32f* pSrcLow, const Ipp32f* pSrcHigh, int srcLen, Ipp16u* pDst, IppsWTInvState_32f16u* pState) {
			return ippsWTInv_32f16u(pSrcLow, pSrcHigh, srcLen, pDst, pState);
		}

		inline auto WTInv(const Ipp32f* pSrcLow, const Ipp32f* pSrcHigh, int srcLen, Ipp32f* pDst, IppsWTInvState_32f* pState) {
			return ippsWTInv_32f(pSrcLow, pSrcHigh, srcLen, pDst, pState);
		}

		inline auto ReplaceNAN_32f(Ipp32f* pSrcDst, int len, Ipp32f value) {
			return ippsReplaceNAN_32f_I(pSrcDst, len, value);
		}

		inline auto ReplaceNAN_64f(Ipp64f* pSrcDst, int len, Ipp64f value) {
			return ippsReplaceNAN_64f_I(pSrcDst, len, value);
		}

		inline auto PatternMatch(const Ipp8u* pSrc, int srcStep, int srcLen, const Ipp8u* pPattern, int patternStep, int patternLen, int patternSize, Ipp16u* pDst, IppPatternMatchMode hint, Ipp8u* pBuffer) {
			return ippsPatternMatch_8u16u(pSrc, srcStep, srcLen, pPattern, patternStep, patternLen, patternSize, pDst, hint, pBuffer);
		}

		inline auto TopKInit(Ipp32s* pDstValue, Ipp64s* pDstIndex, Ipp64s dstLen) {
			return ippsTopKInit_32s(pDstValue, pDstIndex, dstLen);
		}

		inline auto TopKInit(Ipp32f* pDstValue, Ipp64s* pDstIndex, Ipp64s dstLen) {
			return ippsTopKInit_32f(pDstValue, pDstIndex, dstLen);
		}

		inline auto TopK(const Ipp32s* pSrc, Ipp64s srcIndex, Ipp64s srcStride, Ipp64s srcLen, Ipp32s* pDstValue, Ipp64s* pDstIndex, Ipp64s dstLen, IppTopKMode hint, Ipp8u* pBuffer) {
			return ippsTopK_32s(pSrc, srcIndex, srcStride, srcLen, pDstValue, pDstIndex, dstLen, hint, pBuffer);
		}

		inline auto TopK(const Ipp32f* pSrc, Ipp64s srcIndex, Ipp64s srcStride, Ipp64s srcLen, Ipp32f* pDstValue, Ipp64s* pDstIndex, Ipp64s dstLen, IppTopKMode hint, Ipp8u* pBuffer) {
			return ippsTopK_32f(pSrc, srcIndex, srcStride, srcLen, pDstValue, pDstIndex, dstLen, hint, pBuffer);
		}
	}

	inline namespace ipps_l {
		inline auto SortRadixGetBufferSize(IppSizeL len, IppDataType dataType, IppSizeL* pBufferSize) {
			return ippsSortRadixGetBufferSize_L(len, dataType, pBufferSize);
		}

		inline auto SortRadixIndexGetBufferSize(IppSizeL len, IppDataType dataType, IppSizeL* pBufferSize) {
			return ippsSortRadixIndexGetBufferSize_L(len, dataType, pBufferSize);
		}

		inline auto SortRadixAscend_32s_I(Ipp32s* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_32s_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_32f_I(Ipp32f* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_32f_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64u_I(Ipp64u* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64u_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64s_I(Ipp64s* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64s_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixAscend_64f_I(Ipp64f* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixAscend_64f_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_32s_I(Ipp32s* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_32s_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_32f_I(Ipp32f* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_32f_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64u_I(Ipp64u* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64u_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64s_I(Ipp64s* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64s_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixDescend_64f_I(Ipp64f* pSrcDst, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixDescend_64f_I_L(pSrcDst, len, pBuffer);
		}

		inline auto SortRadixIndexAscend_64s(const Ipp64s* pSrc, IppSizeL srcStrideBytes, IppSizeL* pDstIndx, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_64s_L(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexAscend_64u(const Ipp64u* pSrc, IppSizeL srcStrideBytes, IppSizeL* pDstIndx, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexAscend_64u_L(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend_64s(const Ipp64s* pSrc, IppSizeL srcStrideBytes, IppSizeL* pDstIndx, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_64s_L(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}

		inline auto SortRadixIndexDescend_64u(const Ipp64u* pSrc, IppSizeL srcStrideBytes, IppSizeL* pDstIndx, IppSizeL len, Ipp8u* pBuffer) {
			return ippsSortRadixIndexDescend_64u_L(pSrc, srcStrideBytes, pDstIndx, len, pBuffer);
		}
	}
}