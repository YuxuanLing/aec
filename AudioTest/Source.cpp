#include <iostream>
#include <fstream>
#include "pp20aec.h"

#define N 1



void AEC(float *playout_buffer,
	float *capture_buffer,
	float *aec_playout_buffer,
	float *aec_capture_buffer,
	int samples,
	float delay_ms, PP20AEC *aec)
{
	//PP20AEC aec = PP20AEC();
	aec->SetDelay(delay_ms);
	aec->Process(playout_buffer, capture_buffer, aec_playout_buffer, aec_capture_buffer, samples);
}

/*
samples: 480
channels: 1
meta_tag_size: 3
*/
void AECandAGCTest(int samples = 480, int channels = 1, int meta_tag_size = 3)
{
	float audioioDelayMs = 20.0f;
	char ch[6] = { '0', '1', '2', '3', '4' ,'5'};
	PP20AEC aec = PP20AEC();

	for (int k = 1; k <= N; k++)
	{
		std::ifstream myInFromNetwork;
		std::ifstream myInFromMic;
		std::ofstream myOutToNetwork;

		std::string basedir = "E:/personal/coding/aec/13_10_2020_11_09_58 Call Number 1 Audio Dump/13_10_2020_11_09_58 Call Number 1 Audio Dump/13_10_2020_11_09_58";
		//std::string basedir = "E:/personal/coding/aec/13_10_2020_11_17_53 Call number 2 Audio Dump/13_10_2020_11_17_53 Call number 2 Audio Dump";
		std::string my_in_from_network_name = basedir + "/mInFromNetwork.raw";
		std::string my_in_from_mic_name = basedir + "/mInFromMic.raw";

		myInFromNetwork.open(my_in_from_network_name.c_str(), std::ios::binary);
		myInFromMic.open(my_in_from_mic_name.c_str(), std::ios::binary);

		myInFromMic.seekg(0, std::ios::end);
		int lengthMic = myInFromMic.tellg();
		myInFromMic.seekg(0, std::ios::beg);

		int samples_bytes = samples * sizeof(float);
		float *in_buffer_my = (float *)_mm_malloc(samples_bytes, 64);
		float *capture_buffer_my = (float *)_mm_malloc(samples_bytes, 64);
		float *playout_buffer = (float *)_mm_malloc(samples_bytes, 64);
		float *aec_in_buffer = (float *)_mm_malloc(samples_bytes, 64);
		float *out_buffer = (float *)_mm_malloc(samples_bytes, 64);

		//audioioDelayMs = k * 20;
		std::string my_out_to_network_name = (std::string)"./file" + "/myOutToNetwork_xyk" + ch[k] + ".raw";
		myOutToNetwork.open(my_out_to_network_name.c_str(), std::ios::binary);

		// for test
		std::string test_read_mic_name = (std::string)"./file" + "/testReadMic" + ch[k] + ".raw";
		std::ofstream test_read_mic;
		test_read_mic.open(test_read_mic_name.c_str(), std::ios::binary);

		int times = lengthMic / samples_bytes;
		//int times = 1;
		printf("times = %d \n", times);
		for (int i = 0; i < times; i++)
		{
			if (myInFromNetwork.is_open())
			{
				myInFromNetwork.read(reinterpret_cast<char *>(in_buffer_my), samples * sizeof(float));
			}

			if (myInFromMic.is_open())
			{
				myInFromMic.read(reinterpret_cast<char *>(capture_buffer_my), samples * sizeof(float));

				test_read_mic.write(reinterpret_cast<char*>(capture_buffer_my), samples * sizeof(float));
			}

			// Allocate a new aec buffer that gets a mix of all the buffers, and
			// stripped of its meta_tag header (used for 3D)
			//memset(aec_in_buffer, 0, samples_bytes);
			//memcpy(aec_in_buffer, capture_buffer_my, samples_bytes);
			//{
			//	for (int s = 0; s < samples; s++)
			//		for (int b = 0; b < channels; b++)
			//			aec_in_buffer[s] += in_buffer_my[(s + meta_tag_size) * channels + b];
			//}
			//getchar();
			AEC(aec_in_buffer, capture_buffer_my, playout_buffer, out_buffer, samples, audioioDelayMs, &aec);
			//we clip the samples, to make sure they are between -1 and 1
			/* if (Clip(out_buffer, samples) > (0.9 * samples))
			{
			LOG_WARNING_F("Encountered severe clipping, resetting AEC.");
			#if !defined(TARGET_OS_IPHONE) || TARGET_OS_IPHONE == 0
			mAec.Reset();
			#endif
			}*/

			// Audio Gain Control
			// AGC(out_buffer, samples);

			if (myOutToNetwork.is_open())
			{
				myOutToNetwork.write(reinterpret_cast<char *>(out_buffer), samples * sizeof(float));
			}
			printf("times %5d Finished \n",i);
		}

		_mm_free(in_buffer_my);
		_mm_free(capture_buffer_my);
		_mm_free(playout_buffer);
		_mm_free(aec_in_buffer);
		_mm_free(out_buffer);

		myInFromMic.close();
		myInFromNetwork.close();
		myOutToNetwork.close();
	}
	printf("test Finished \n");
}

int main()
{
	AECandAGCTest();
	getchar();
	return 0;
}