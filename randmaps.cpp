#include "globalvars.h"

data_t func1(data_t x)
{
	IppStatus status;
	int srclen;
	srclen = 16;
	Ipp8u pSrc[16];
	Ipp8u pDst[16];

	x &= mask;
	data_t y;
	for (size_t i = 0; i < 16; i++)
	{
		pSrc[i] = 0;
		pDst[i] = 0;
	}
#if defined(_MSC_VER)
	memcpy(pSrc, &x, sizeof(data_t));
#else
	for (size_t i = 0; i < sizeof(data_t); i++)
	{
		pSrc[i] = ((u8*)(&x))[i];
	}
#endif
	status = ippsAESEncryptECB(pSrc, pDst, srclen, pCtx1);
	switch (status)
	{
	case ippStsNoErr: /*cout << "ippStsNoErr: Indicates no error." << endl;*/ break;
	case ippStsNullPtrErr: cout << "ippStsNullPtrErr: Indicates an error condition if the specified  pointer is NULL." << endl; break;
	case ippStsLengthErr: cout << "ippStsLengthErr: Indicates an error condition if the input data stream length is less than or equal to zero." << endl; break;
	case ippStsUnderRunErr: cout << "ippStsUnderRunErr: Indicates an error condition if srclen is not divisible by cipher block size." << endl; break;
	case ippStsContextMatchErr: cout << "ippStsContextMatchErr: Indicates an error condition if the context parameter does not match the operation." << endl; break;
	default: break;
	}
#if defined(_MSC_VER)
	memcpy(&y, pDst, sizeof(data_t));
#else
	for (size_t i = 0; i < sizeof(data_t); i++)
	{
		((u8*)(&y))[i] = pDst[i];
	}
#endif
	return (y & mask);
}

data_t func2(data_t x)
{
	IppStatus status;
	Ipp8u pSrc[16];
	Ipp8u pDst[16];
	int srclen;
	srclen = 16;

	x &= mask;
	data_t y;
	for (size_t i = 0; i < 16; i++)
	{
		pSrc[i] = 0;
		pDst[i] = 0;
	}
#if defined(_MSC_VER)
	memcpy(pSrc, &x, sizeof(data_t));
#else
	for (size_t i = 0; i < sizeof(data_t); i++)
	{
		pSrc[i] = ((u8*)(&x))[i];
	}
#endif
	status = ippsSMS4DecryptECB(pSrc, pDst, srclen, pCtx2);
	switch (status)
	{
	case ippStsNoErr: /*cout << "ippStsNoErr: Indicates no error." << endl;*/ break;
	case ippStsNullPtrErr: cout << "ippStsNullPtrErr: Indicates an error condition if the specified  pointer is NULL." << endl; break;
	case ippStsLengthErr: cout << "ippStsLengthErr: Indicates an error condition if the input data stream length is less than or equal to zero." << endl; break;
	case ippStsUnderRunErr: cout << "ippStsUnderRunErr: Indicates an error condition if srclen is not divisible by cipher block size." << endl; break;
	case ippStsContextMatchErr: cout << "ippStsContextMatchErr: Indicates an error condition if the context parameter does not match the operation." << endl; break;
	default: break;
	}
#if defined(_MSC_VER)
	memcpy(&y, pDst, sizeof(data_t));
#else
	for (size_t i = 0; i < sizeof(data_t); i++)
	{
		((u8*)(&y))[i] = pDst[i];
	}
#endif
	return (y & mask);
}