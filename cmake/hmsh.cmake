if(TARGET hmsh)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
        hmsh
        GIT_REPOSITORY https://github.com/komietty/hmsh
        GIT_TAG main
)
FetchContent_MakeAvailable(hmsh)

