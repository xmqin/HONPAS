! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module jobList

  implicit none

PUBLIC:: &
  countJobs,  &! count number of jobs and lists
  runJobs,    &! runs a job list read from a file unit
  getResults   ! collect requested results into summary files

PRIVATE    ! nothing is declared public beyond this point

  ! Internal parameters
  character(len=*),parameter:: dataSuffix = '.fdf'
  character(len=*),parameter:: endSuffix = '.EIG'
  character(len=*),parameter:: file0Exit = '0_NORMAL_EXIT'
  character(len=*),parameter:: defaultFiles = &
                                       '*.fdf *.vps *.psf *.ion queue.sh'
  character(len=*),parameter:: defaultQueue = &
                                       './siesta < $jobName.fdf > $jobName.out'
  character(len=*),parameter:: defaultRequest(1) = (/'energy'/)
  integer,parameter:: maxLines = 1000   ! Max lines in job-list file
  integer,parameter:: maxWords = 100    ! Max words in one line
  integer,parameter:: maxJobs = 1000    ! Max jobs in job list
  integer,parameter:: ll      = 500     ! Max characters per line
  integer,parameter:: wl      = 500     ! Max characters per word
  integer,parameter:: unitIn = 42       ! I/O unit for input files
  integer,parameter:: unitOut = 43      ! I/O unit for output files
  integer,parameter:: dp = kind(1.d0)

  ! Internal variables and arrays
  integer,save:: jobCores, totCores, totJobs, totLists

  ! Derived type to hold calculation results
  type resultsType
    private
    real(dp)        :: cell(3,3)   =0._dp
    real(dp)        :: stress(3,3) =0._dp
    real(dp)        :: volume      =0._dp
    real(dp)        :: energy      =0._dp
    real(dp)        :: pressure    =0._dp
    real(dp)        :: virial      =0._dp 
    real(dp),pointer:: coord(:,:)  =>null()
    real(dp),pointer:: force(:,:)  =>null()
    integer, pointer:: za(:)       =>null()   ! atomic numbers
  end type resultsType

CONTAINS

!------------------------------------------------------------------------------

subroutine countJobs( unit, nLists, nJobs, nCores )

! Counts number of jobs and lists in a joblist file, specified by a I/O unit

  implicit none
  integer,         intent(in) :: unit   ! IO unit of datafile
  integer,optional,intent(out):: nLists ! total number of job lists
  integer,optional,intent(out):: nJobs  ! total number of jobs
  integer,optional,intent(out):: nCores ! total number of required cores

  integer :: nc
  character(len=wl):: myDir

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(unit,myDir,defaultQueue,defaultFiles,defaultRequest,' ','count')

  ! Copy number of jobs and lists to output variables
  if (present(nLists)) nLists = totLists
  if (present(nJobs))  nJobs  = totJobs
  if (present(nCores)) nCores = totCores

end subroutine countJobs

!------------------------------------------------------------------------------

subroutine runJobs( unit )

! Runs or queues a number of jobs in a joblist file, specified by a I/O unit

  implicit none
  integer,intent(in) :: unit   ! IO unit of datafile

  integer :: nc
  character(len=wl):: myDir

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(unit,myDir,defaultQueue,defaultFiles,defaultRequest,' ','run')

end subroutine runJobs

!------------------------------------------------------------------------------

subroutine getResults( unit )

! Collects requested results of jobs in a joblist file (specified by a 
! I/O unit), into summary output files, within each sub-list directory

  implicit none
  integer,intent(in) :: unit   ! IO unit of datafile

  integer :: nc
  character(len=wl):: myDir

  ! Get current directory and add trailing '/' if necessary
  call getcwd(myDir)
  nc = len(trim(myDir))
  if (myDir(nc:nc) /= '/') myDir = myDir(1:nc) // '/'

  ! Scan the 'root' list (the whole file), specified by blank name
  call scanList(unit,myDir,defaultQueue,defaultFiles,defaultRequest,' ','get')

end subroutine getResults

!------------------------------------------------------------------------------

recursive subroutine scanList( unit, dir, queue, files, request, &
                               listName, task )

  implicit none
  integer,         intent(in) :: unit       ! I/O unit of joblist file
  character(len=*),intent(in) :: dir        ! parent directory
  character(len=*),intent(in) :: queue      ! queuing statement
  character(len=*),intent(in) :: files      ! required data files
  character(len=*),intent(in) :: request(:) ! requested results
  character(len=*),intent(in) :: listName   ! name of job list
  character(len=*),intent(in) :: task       ! ('count'|'run'|'get')

  character(len=1),parameter:: separator(1) = (/' '/)
  integer :: iCase, iLine, iostat, nCases, nCores, nResults, nWords, nJobs
  character(len=wl):: caseDir(maxJobs), caseName(maxJobs), caseType(maxJobs), &
                       fileIn, fileOut, line, myDir, myFiles, myQueue, &
                       myRequest(maxWords), newList, words(maxWords)
  real(dp):: results(maxWords)
  logical :: finished

  ! Set the list directory and copy data files to it
  if (listName==' ') then  ! this is the 'root' list (the whole file)
    myDir = dir
  else                     ! this is a genuine list of jobs
    myDir = trim(dir) // trim(listName) // '/'
    ! Create a new directory for this list and copy files to it
    if (task=='run') then
      call system('mkdir -p ' // trim(myDir))
      call system('cp -f ' // trim(files) // ' ' // trim(myDir) // &
                  ' 2> /dev/null' )
      call chdir(trim(myDir))
    endif
  endif ! (listName==' ')

  ! Global initializations
  if (listName==' ') then   ! this is the 'root' list (the whole file)
    totCores = 0            ! total number of cores
    totJobs = 0             ! total number of jobs
    totLists = 0            ! total number of job lists
  endif

  ! Set my own copies of queue, files, and request
  nResults = count(request/=' ')
  myRequest = ' '
  myRequest(1:nResults) = request(1:nResults)
  myQueue = queue
  myFiles = files
  call getCores(myQueue,jobCores)   ! get number of cores per job

  ! Loop on lines of datafile
  finished = .false.       ! has the list terminated normally?
  nJobs = 0                ! number of jobs in list
  nCores = 0               ! number of cores required by jobs in list
  nCases = 0               ! number of 'cases' (jobs or sublists) in list
  do iLine = 1,maxLines

    ! Read one line (possibly continued) of input file, ignoring comments
    call readLine( unit, line, iostat )

    ! End of file check
    if (iostat<0) then                      ! end of file
      if (listName==' ') then               ! this is the 'root' list
        finished = .true.
        exit ! iLine loop
      else
        print*,'runList ERROR: unterminated %list = '//trim(listName)
        stop
      endif ! (listName==' ')
    endif ! (iostat<0)

    ! Parse line
    call parser(separator,line,words,nWords)

    ! Act depending on line content
    if (nWords==0) then                     ! blank or comment line
      cycle
    elseif (trim(words(1))=='%queue') then  ! queue specification
      line = adjustl(line)
      myQueue = adjustl(line(7:))           ! remove '%queue' from line
      call getCores(myQueue,jobCores)
    elseif (trim(words(1))=='%files') then  ! data files specification
      line = adjustl(line)
      myFiles = adjustl(line(7:))           ! remove '%files' from line
    elseif (trim(words(1))=='%result') then ! result-request especification
      myRequest = ' '
      myRequest(1:nWords-1) = words(2:nWords) ! exclude 1st word (it is %result)
    elseif (words(1)=='%list') then         ! new sub-list
      nCases = nCases+1
      newList = words(2)
      caseType(nCases) = 'list'
      caseName(nCases) = newList
      caseDir(nCases) = trim(myDir) // trim(newList) // '/'
      call scanList(unit,myDir,myQueue,myFiles,myRequest,newList,task)
    elseif (words(1)=='%endlist') then      ! end of my list
      if (words(2)==listName) then
        finished = .true.
        exit ! iLine loop
      else
        print*,'runList ERROR: inconsistent %list = '//trim(listName)// &
               ', %endList = '//trim(words(2))
        stop
      endif
    else                                    ! new job
      nCases = nCases+1
      nJobs = nJobs+1
      nCores = nCores + jobCores
      caseType(nCases) = 'job'
      call nameJob( line, caseName(nCases) )
      caseDir(nCases) = trim(myDir) // trim(caseName(nCases)) // '/'
      if (task=='run') call runOneJob(myDir,myQueue,myFiles,line)
    endif

  end do ! iLine
  if (.not.finished) stop 'runList ERROR: parameter maxLines too small'

  ! Write file of requested results
  if (task=='get') then
    if (listName==' ') then
      fileOut = trim(myDir)//'jobList.results'
    else
      fileOut = trim(myDir)//trim(listName)//'.results'
    endif
    call system('rm -f '//trim(fileOut))  ! clear output file, if it exists
    do iCase = 1,nCases
      if (caseType(iCase)=='list') then
        fileIn = trim(caseDir(iCase))//trim(caseName(iCase))//'.results'
        call system('cat '//trim(fileIn)//' >> '//trim(fileOut))
        call system('echo '' '' >> '//trim(fileOut))
      else
        call readResult( caseDir(iCase), myRequest, caseName(iCase), &
                         results=results )
        open(unitOut,file=trim(fileOut),position='append',status='unknown')
        write(unitOut,'(10e18.9)') results(1:nResults) 
        close(unitOut)
      endif
    end do
  endif ! (task=='get')

  ! Count jobs and list (but not if it is a list of lists)
  totJobs = totJobs+nJobs
  totCores = totCores+nCores
  if (nJobs>0) totLists = totLists+1

end subroutine scanList

!------------------------------------------------------------------------------

subroutine getCores( queue, ncores )

! Extracts an integer number from string 'queue', assuming that it is the
! number of cores per job, or returns 1 if there is no such integer

  implicit none
  character(len=*),intent(in) :: queue    ! queue statement string
  integer,         intent(out):: ncores   ! number of cores per job

  character(len=1),parameter:: separator(1) = (/' '/)
  character(len=wl):: words(maxWords)
  integer:: iWord, n, nWords

  ! Parse queue into blank-separated words
  call parser(separator,queue,words,nWords)

  ! Find a word that is an integer number. Assume that it is the num. of cores
  ncores = 1                               ! default value
  do iWord = 1,nWords
    n = verify(trim(words(iWord)),'0123456789')  ! verify is an intrinsic funct
    if (n==0) then                         ! all characters are numbers
      read(words(iWord),*) ncores          ! read ncores from word
      return
    endif
  end do

end subroutine getCores

!------------------------------------------------------------------------------

subroutine nameJob( jobLine, jobName )

! Generates a job name, concatenating words in a job-specification line
! .fdf suffixes are removed from words

  implicit none
  character(len=*),intent(in) :: jobLine   ! string of job specifications
  character(len=*),intent(out):: jobName   ! job name

  character(len=1),parameter:: separators(2) = (/' ',';'/)
  integer:: iWord, nc, nWords
  character(len=wl):: line, word, words(maxWords)

  ! Parse job line
  line = jobLine
  call parser(separators,line,words,nWords)

  ! Set job name by concatenating words
  jobName = ' '
  do iWord = 1,nWords
    word = words(iWord)
    nc = len(trim(word))                         ! number of characters in word
    if (word(nc-3:nc)=='.fdf') word=word(1:nc-4) ! remove '.fdf' from word
    jobName = trim(jobName)//'_'//word           ! add word to job name
  end do ! iWord
  jobName=jobName(2:)                            ! remove leading '_'

end subroutine nameJob

!------------------------------------------------------------------------------

subroutine runOneJob( dir, queue, files, jobLine )

! Submits (queues) one job, specified by string 'jobLine' of specifications,
! using the 'queue' statement, within a job-specific subdirectory of 'dir'
! The subdirectory is named by concatenating the words in 'jobLine' (without
! '.fdf' suffixes), separated by '_'

  implicit none
  character(len=*),intent(in) :: dir       ! work directory
  character(len=*),intent(in) :: queue     ! queuing statement
  character(len=*),intent(in) :: files     ! data files
  character(len=*),intent(in) :: jobLine   ! line of job specification

  character(len=1),parameter:: separator(1) = (/';'/)
  logical:: jobEnded
  integer:: ic, iLine, iWord, jc, jWord, nc, nLines, nWords
  character(len=wl):: fileName, jobDir, jobName, line(maxWords), &
                      myLine, queueJob, word, words(maxWords)

  ! Find job name
  call nameJob( jobLine, jobName )

  ! Create a new directory for this job and copy files to it
  jobDir = trim(dir) // trim(jobName) // '/'
  call system('mkdir -p ' // trim(jobDir))
  call system('cp -f ' // trim(files) // ' ' // trim(jobDir) // &
              ' 2> /dev/null')

  ! Parse job-specifications line
  myLine = jobLine
  call parser(separator,myLine,words,nWords)

  ! Convert job specifications to fdf format lines
  nLines = nWords                      ! lines in the output fdf file
  do iWord = 1,maxWords                ! loop on job specifications
    word = words(iWord)
    nc = len(trim(word))
    if (word(nc-3:nc)=='.fdf') then    ! include new .fdf file
      line(iWord) = '%include '//trim(word)
    else                               ! copy statement to new line
      line(iWord) = word
    end if
  end do ! iWord

  ! Write .fdf file for job, but avoid overwriting it (this would occur only
  ! if job contains a single .fdf file, with no other specifications)
  fileName = trim(jobName)//'.fdf'
  if (trim(fileName)/=trim(words(1))) then
    open(unitOut,file=trim(jobDir)//trim(fileName))
    write(unitOut,'(a)') 'SystemLabel '//trim(jobName) ! this is siesta-specific
    do iLine = nLines,1,-1               ! last lines (words) have priority
      write(unitOut,'(a)') trim(line(iLine))
    end do
    close(unitOut)
  endif

  ! Set job
  queueJob = queue
  nc = len(trim(queueJob))
  jc = len('$jobName')
  do ic = nc-jc+1,1,-1
    if (queueJob(ic:ic+jc-1)=='$jobName') &
      queueJob(ic:) = trim(jobName)//queueJob(ic+jc:)
  end do
!  print*, 'runOneJob: queueJob = ',trim(queueJob)

  ! Submit job only if directory does not contain terminated results already
  ! This is intended to re-run a whole list for failed jobs

  call chdir(trim(jobDir))

  ! Check for existence of a special file which signals the end of execution
  inquire(file=trim(file0Exit),exist=jobEnded)

  if (.not.jobEnded) then
     ! Fall back to old test with EIG file, in case of an older
     ! version of Siesta without the new flag file.
     inquire(file=trim(jobName)//endSuffix,exist=jobEnded)
  endif

  if (.not.jobEnded) then
     !print "(2a)", 'queueJob = ',trim(queueJob)
     call system(trim(queueJob))
  endif
  call chdir(trim(dir))

end subroutine runOneJob

!------------------------------------------------------------------------------

subroutine readLine( unit, line, iostat )

! Reads one line of a job list file (specified by a I/O unit), removing 
! comments (marked by a leading '#'), and including continuation lines 
! (marked by a trailing '\')

  implicit none
  integer,         intent(in) :: unit    ! input file unit
  character(len=*),intent(out):: line    ! line read
  integer,optional,intent(out):: iostat  ! I/O status

  logical keepReading
  integer:: lc, nc, newc, status
  character:: lastChar
  character(len=1024):: newLine    ! a sufficiently long array

  line = ' '                       ! set a blank line to start
  nc = 0                           ! number of characters in line
  keepReading = .true.
  do while (keepReading)           ! loop on continuation lines

    ! Read new line from input file
    read(unit,'(a)',iostat=status) newLine

    ! End-of-file trap
    if (present(iostat)) iostat=status
    if (status<0) then
      line = ' '
      return
    endif

    ! Remove comments and leading blanks
    lc = scan(newLine,'#')-1             ! last character before comments
    if (lc>=0) newLine = newLine(1:lc)   ! select characters before comments
    newLine = adjustl(newLine)           ! remove leading blanks

    ! Find last character and remove continuation mark
    newc = len(trim(newLine))            ! nonblank characters in new line
    lastChar = newLine(newc:newc)        ! last nonblank character
    if (lastChar==char(92)) then         ! line will continue  (backslash) '\'
      newLine(newc:newc) = ' '           ! remove '\' mark
      newc = len(trim(newLine))          ! remaining nonblank characters
      keepReading = .true.
    else
      keepReading = .false.
    endif

    ! Concatenate continued line
    if (nc==0) then                      ! this is first line
      line = newLine(1:newc)
      nc = newc
    else                                 ! this is a continuation line
      line = line(1:nc)//' '//newLine(1:newc)
      nc = nc+1+newc
    endif

    ! Check array length
    if (nc>len(line)) stop 'readLine ERROR: len(line) too small'

  end do ! while (keepReading)

end subroutine readLine

!------------------------------------------------------------------------------

subroutine parser( separators, line, words, nWords )

! Parses one line into words, using any of a number of separators.
! If ' ' is not one of the separators, the 'words' may contain blanks,
! but leading and trailing blanks will be removed from them in any case

  implicit none
  character(len=*),intent(in)   :: separators(:) ! character(s) between words
  character(len=*),intent(in)   :: line      ! line to be parsed
  character(len=*),intent(out)  :: words(:)  ! words in line, without comments
  integer,         intent(out)  :: nWords    ! number of words

  character(len=ll):: myLine
  integer:: is, iWord, nc, ncs, ns

  ! Loop on words
  myLine = adjustl(line)                  ! copy line, removing leading blanks
  ns = size(separators)                   ! number of alternative separators
  do iWord = 1,size(words)                ! avoid overflooding array size
    nc = len(trim(myLine))                ! number of nonblank characters
    do is = 1,ns                          ! loop on separators
      ncs = scan(myLine,separators(is))-1 ! last character before this separator
      if (ncs>0) nc = min(nc,ncs)         ! character before first separator
    end do
    if (nc<=0) then                       ! no more words
      nWords = iWord-1  
      return                              ! normal return point
    else
      if (nc>len(words(iWord))) stop 'parser ERROR: len(words) too small'
      words(iWord) = adjustl(myLine(1:nc))
      myLine = adjustl(myLine(nc+2:))     ! discard previous word from myLine
    endif
  end do ! iWord
  stop 'parser ERROR: size(words) too small'

end subroutine parser

!------------------------------------------------------------------------------

subroutine readResult( dir, request, name, result, results )

! Reads magnitudes specified by 'request' string(s), calculated by a SIESTA 
! job named 'name', within directory 'dir'. The magnitudes are output in one
! or both of: derived-type structure 'result'; and real array 'results'
! Presently allowed values (case sensitive) of 'request' are
!   (volume|energy|pressure|virial|maxForce|avgForce)
! Additionally, derived-type 'result' returns the system geometry, forces,
! and stress tensor

  implicit none
  character(len=*),          intent(in) :: dir         ! job directory
  character(len=*),          intent(in) :: request(:)  ! requested resuts
  character(len=*),          intent(in) :: name        ! job name
  type(resultsType),optional,intent(out):: result      ! job-results structure
  real(dp),         optional,intent(out):: results(:)  ! requested job results

  integer :: ia, ic, is, iostat, iResult, nAtoms, nResults, za
  real(dp):: c(3,3)
  character(len=wl):: fileName
  type(resultsType) :: res

  ! Read final structure
  fileName = trim(dir)//trim(name)//'.XV'
  open(unitIn,file=trim(fileName),status='old',iostat=iostat)
!  print'(a,i3,2x,a)','readResult: iostat,fileName=', iostat, trim(fileName)
  if (iostat==0) then
    do ic = 1,3
      read(unitIn,*) res%cell(:,ic)
    end do
    read(unitIn,*) nAtoms
    allocate( res%za(nAtoms), res%coord(3,nAtoms) )
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%coord(:,ia)
    end do
  else
    nAtoms = 0
    allocate( res%za(nAtoms), res%coord(3,nAtoms) )
    res%cell = 0
    res%coord = 0
    res%za = 0
  end if
  close(unitIn)
  
  ! Read energy, force, and stress
  fileName = trim(dir)//'FORCE_STRESS'
  open(unitIn,file=trim(fileName),status='old',iostat=iostat)
  if (iostat==0) then
    read(unitIn,*) res%energy
    read(unitIn,*) res%stress
    read(unitIn,*) nAtoms
    deallocate(res%za)
    allocate(res%force(3,nAtoms),res%za(nAtoms))
    do ia = 1,nAtoms
      read(unitIn,*) is, res%za(ia), res%force(:,ia)
    end do
  else
    allocate( res%force(3,0) )
    res%energy = 0
    res%stress = 0
    res%force = 0
  end if
  close(unitIn)

  ! Find cell volume, pressure, and virial
  c = res%cell
  res%volume = abs( c(1,1)*c(2,2)*c(3,3) - c(1,1)*c(2,3)*c(3,2) + &
                    c(1,2)*c(2,3)*c(3,1) - c(1,2)*c(2,1)*c(3,3) + &
                    c(1,3)*c(2,1)*c(3,2) - c(1,3)*c(2,2)*c(3,1) )
  res%pressure = -res%stress(1,1)-res%stress(2,2)-res%stress(3,3)
  res%virial = sum(res%coord*res%force)

  ! Copy results to output structure
  if (present(result)) then
    result = res
  end if

  ! Select requested (scalar) results
  if (present(results)) then
    nResults = count(request/=' ')
    do iResult = 1,nResults
      select case (trim(request(iResult)))
        case ('energy')
          results(iResult) = res%energy
        case ('volume')
          results(iResult) = res%volume
        case ('pressure')
          results(iResult) = res%pressure
        case ('virial')
          results(iResult) = res%virial
        case ('maxForce')
          results(iResult) = sqrt(maxval(sum(res%force**2,1)))
        case ('avgForce')
          results(iResult) = sqrt(sum(res%force**2)/max(1,nAtoms))
      end select
    end do
  end if

  ! Deallocate internal pointers
  if (associated(res%za)) deallocate(res%za)
  if (associated(res%coord)) deallocate(res%coord)
  if (associated(res%force)) deallocate(res%force)

end subroutine readResult

END MODULE jobList
