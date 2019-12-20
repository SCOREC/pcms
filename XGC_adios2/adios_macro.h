#ifdef ADIOS
#define ADIOS_OPEN(fd,grp,file,mode,grp_comm,err)\
  call adios_open(fd,trim(grp)//char(0), \
       trim(file)//char(0), \
       mode//char(0),grp_comm,err)
#define ADIOS_SET_PATH(fd,path,err)			\
  call adios_set_path(fd,trim(path)//char(0),err)
#define ADIOS_SET_PATH_VAR(fd,path,var,err)				\
  call adios_set_path_var(fd,trim(path)//char(0),trim(var)//char(0),err)
#define ADIOS_GROUP_SIZE(fd,size,total_size,err)	\
  call adios_group_size(fd,size,total_size,err)

#if defined(__PGI) && (__PGIC__<15) || defined(__GFORTRAN__)
#define ADIOS_WRITE(fd,var,err)			\
  call adios_write(fd,'var'//char(0),var,err)
#else
#define ADIOS_WRITE(fd,var,err)			\
  call adios_write(fd,#var//char(0),var,err)
#endif

#define ADIOS_WRITE_LBL(fd,label,var,err)		\
  call adios_write(fd,trim(label)//char(0),var,err)

#if defined(__PGI) && (__PGIC__<15) || defined(__GFORTRAN__)
#define ADIOS_READ(fd,var,size,err)			\
  call adios_read(fd,'var'//char(0),var,size,err)
#else
#define ADIOS_READ(fd,var,size,err)			\
  call adios_read(fd,#var//char(0),var,size,err)
#endif
#define ADIOS_READ_LBL(fd,label,var,size,err)		\
  call adios_read(fd,trim(label)//char(0),var,size,err)
#define ADIOS_CLOSE(fd,err)			\
  call adios_close(fd,err)
#endif
