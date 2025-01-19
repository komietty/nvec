if(TARGET hmesh)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
        hmsh
        GIT_REPOSITORY https://github.com/komietty/hmesh
        GIT_TAG main
)
FetchContent_MakeAvailable(hmsh)

