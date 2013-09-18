#ifndef TIMER_HPP
#define TIMER_HPP



#if defined WIN32

    #include <windows.h>

#elif defined __APPLE__

    #include <sys/time.h>
    #include <mach/mach_time.h> 

#else //Linux

    #include <time.h>

#endif



namespace zipHMM {

    class Timer {

#if defined WIN32
  
        LARGE_INTEGER frequency;
        LARGE_INTEGER startTime, endTime;

    public:

        Timer()
        {
            QueryPerformanceFrequency(&frequency);
        } 
   
        void start()
        {
            QueryPerformanceCounter(&startTime);
        }
    
        void stop()
        {
            QueryPerformanceCounter(&endTime);
        }
    
        double timeElapsed() const
        {
            return double(endTime.QuadPart - startTime.QuadPart) / frequency.QuadPart;
        }
    
#elif defined __APPLE__

        uint64_t startTime, endTime;
        mach_timebase_info_data_t info;

    public :
        Timer() {
            mach_timebase_info(&info);
        }

        void start()
        {
            startTime = mach_absolute_time();
            startTime *= info.numer;
            startTime /= info.denom;
        }

        void stop()
        {
            endTime = mach_absolute_time();
            endTime *= info.numer;
            endTime /= info.denom;
        }


        double timeElapsed() const
        {
            return (endTime - startTime) / 1000000000.0;
        }
    
#else //Linux

        timespec startTime, endTime;

    public:

        Timer() {} 
   
        void start()
        {
            clock_gettime(CLOCK_REALTIME, &startTime);
        }
    
        void stop()
        {
            clock_gettime(CLOCK_REALTIME, &endTime);
        }
    
        double timeElapsed() const
        {
            timespec diff;
            diff.tv_sec = endTime.tv_sec - startTime.tv_sec;
            diff.tv_nsec = endTime.tv_nsec - startTime.tv_nsec;

            return double(diff.tv_sec) + double(diff.tv_nsec) / 1000000000.0;
        }
    
#endif
    };
}

#endif
