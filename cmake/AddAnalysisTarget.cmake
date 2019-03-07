function(ADD_ANALYSIS target_name)
        message(STATUS "-----------")
        #+++++++++++++++++++++++++
        # CUDA TARGETS           |
        #+++++++++++++++++++++++++
        if(BUILD_CUDA_TARGETS)
                 message(STATUS "Adding target ${target_name} to CUDA backend. Executable file name: ${PROJECT_NAME}-${target_name}_cuda")
                  
                 cuda_add_executable("${target_name}_cuda"
                 #EXCLUDE_FROM_ALL 
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cu"    
                 OPTIONS -Xcompiler -DHYDRA_DEVICE_SYSTEM=CUDA -DHYDRA_HOST_SYSTEM=CPP)
                 
                 set_target_properties( "${target_name}_cuda" 
                 PROPERTIES OUTPUT_NAME "${PROJECT_NAME}-${target_name}_cuda")
                 
                 target_link_libraries("${target_name}_cuda" ${ROOT_LIBRARIES} ${LIBCONFIGPP_LIBRARIES})
                 
                 install(TARGETS "${target_name}_cuda" DESTINATION bin)
          
                
        endif(BUILD_CUDA_TARGETS)
    
        #+++++++++++++++++++++++++
        # TBB TARGETS            |
        #+++++++++++++++++++++++++
        if(BUILD_TBB_TARGETS)
                 message(STATUS "Adding target ${target_name} to TBB backend. Executable file name: ${PROJECT_NAME}-${target_name}_tbb")
                 add_executable("${target_name}_tbb"
                 # EXCLUDE_FROM_ALL
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cpp" )
                    
                 set_target_properties( "${target_name}_tbb" 
                 PROPERTIES COMPILE_FLAGS "-DHYDRA_HOST_SYSTEM=CPP -DHYDRA_DEVICE_SYSTEM=TBB" OUTPUT_NAME "${PROJECT_NAME}-${target_name}_tbb")
                    
                 target_link_libraries( "${target_name}_tbb" ${ROOT_LIBRARIES} ${TBB_LIBRARIES} ${LIBCONFIGPP_LIBRARIES}  )
                 
                 install(TARGETS "${target_name}_tbb" DESTINATION bin)
                   
                       
         endif(BUILD_TBB_TARGETS)
         
        #+++++++++++++++++++++++++
        # CPP TARGETS            |
        #+++++++++++++++++++++++++
        if(BUILD_CPP_TARGETS)
                 message(STATUS "Adding target ${target_name} to CPP backend. Executable file name: ${PROJECT_NAME}-${target_name}_cpp")
                 add_executable("${target_name}_cpp"
                 # EXCLUDE_FROM_ALL 
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cpp" )

                 set_target_properties( "${target_name}_cpp" 
                 PROPERTIES COMPILE_FLAGS "-DHYDRA_HOST_SYSTEM=CPP -DHYDRA_DEVICE_SYSTEM=CPP" OUTPUT_NAME "${PROJECT_NAME}-${target_name}_cpp")
                    
                 target_link_libraries( "${target_name}_cpp" ${ROOT_LIBRARIES} ${TBB_LIBRARIES}  ${LIBCONFIGPP_LIBRARIES} )
                 
                 install(TARGETS "${target_name}_cpp" DESTINATION bin)
              
                       
         endif(BUILD_CPP_TARGETS)
         
          
        #+++++++++++++++++++++++++
        # OMP TARGETS            |
        #+++++++++++++++++++++++++
        if(BUILD_OMP_TARGETS)
                 message(STATUS "Adding target ${target_name} to OMP backend. Executable file name: ${PROJECT_NAME}-${target_name}_omp")
                 add_executable("${target_name}_omp" 
                 #EXCLUDE_FROM_ALL
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cpp" )
                    
                 set_target_properties( "${target_name}_omp" 
                 PROPERTIES COMPILE_FLAGS "-DHYDRA_HOST_SYSTEM=CPP -DHYDRA_DEVICE_SYSTEM=OMP ${OpenMP_CXX_FLAGS}" OUTPUT_NAME "${PROJECT_NAME}-${target_name}_omp")
                    
                 target_link_libraries( "${target_name}_omp" ${ROOT_LIBRARIES} ${OpenMP_CXX_LIBRARIES}  ${LIBCONFIGPP_LIBRARIES} )
                 
                 install(TARGETS "${target_name}_omp" DESTINATION bin)
                   
                       
         endif(BUILD_OMP_TARGETS)
         
endfunction(ADD_ANALYSIS)  
