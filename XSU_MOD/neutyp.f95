! A fortran95 program for G95
! By WQY
MODULE UTILITARIOS_FLOAT
      INTERFACE DIRAPL
        MODULE PROCEDURE DIRACOI
        MODULE PROCEDURE DIRACOF
      END INTERFACE DIRAPL

      CONTAINS
      SUBROUTINE DIRACOF(NWF,IO,DATA,LENGTH,INDEX)
!
!DIRECT ACCESS TO FILE NWF
!IO    = 1/2 WRITE/READ OPTION
!DATA()= MEMORY AREA TO BE WRITTEN OR READ
!LENGTH= AMOUNT OF SINGLE-WORD DATA TO BE TRANSFERRED
!INDEX = START INDEX OF THE !ORRESPONDING AREA ON NWF
!
!USE STANDARD FORTRAN VII STATEMENT WRITE AND READ(NWF,REC=N)
!RECORD LENGTH 128 WORDS (512 BYTES)
!
      INTEGER,INTENT(IN) :: NWF
      INTEGER,INTENT(IN) :: IO
      REAL,DIMENSION(:),INTENT(INOUT) :: DATA
      INTEGER,INTENT(IN) :: LENGTH
      INTEGER,INTENT(INOUT) :: INDEX
      LOGICAL               SWPLUS
      REAL               RECORD(128)
!
      DATA       LREC  / 128 /
!
      SWPLUS=.FALSE.
!    GO TO 10
!   ======================================
!      ENTRY DIRAPL(NWF,IO,DATA,LENGTH,INDEX)
!   ======================================
      SWPLUS=.TRUE.
! ----------------------- INITIATE
   10 NREC =INDEX/LREC
      N     =INDEX-LREC*NREC
      NREC =NREC+1
      LMIN  =1
      LMAX  =LREC-N
      IOUSED=IOUSED+LENGTH
!      GO TO 100
     GO TO (40,100),IO
!
!WRITE. (READ FIRST /LAST RECORD)    ------------------- WRITE(IO=1)
!
  40 IF(N.EQ.0) GO TO 50
     LMAX  =MIN0(LMAX,LENGTH)
     GO TO 60
!!
  50 IF(LMAX.LE.LENGTH) GO TO 70
     LMAX  =LENGTH
  60 READ(NWF,REC=NREC,ERR=70) RECORD
!
  70 DO 80 L=LMIN,LMAX
     N     =N+1
  80 RECORD(N)=DATA(L)
!---------------------------
     WRITE(NWF,REC=NREC) RECORD
!---------------------------
     IF(LMAX.EQ.LENGTH) GO TO 140
     N     =0
     NREC =NREC+1
     LMIN  =LMAX+1
     LMAX  =LMAX+LREC
     GO TO 50
!
!READ-------------------- READ (IO=2)
!
  100 LMAX  =MIN0(LMAX,LENGTH)
!--------------------------
  110 READ(NWF,REC=NREC) RECORD
!--------------------------
      DO 120 L=LMIN,LMAX
      N     =N+1
  120 DATA(L)=RECORD(N)
!
      IF(LMAX.EQ.LENGTH) GO TO 140
      N     =0
      NREC =NREC+1
      LMIN  =LMAX+1
      LMAX  =MIN0(LMAX+LREC,LENGTH)
      GO TO 110
!
  140 IF(SWPLUS) INDEX=INDEX+LENGTH
      RETURN
      END SUBROUTINE DIRACOF

      SUBROUTINE DIRACOI(NWF,IO,DATA,LENGTH,INDEX)
!
!DIRECT ACCESS TO FILE NWF
!IO    = 1/2 WRITE/READ OPTION
!DATA()= MEMORY AREA TO BE WRITTEN OR READ
!LENGTH= AMOUNT OF SINGLE-WORD DATA TO BE TRANSFERRED
!INDEX = START INDEX OF THE !ORRESPONDING AREA ON NWF
!
!USE STANDARD FORTRAN VII STATEMENT WRITE AND READ(NWF,REC=N)
!RECORD LENGTH 128 WORDS (512 BYTES)
!
      INTEGER,INTENT(IN) :: NWF
      INTEGER,INTENT(IN) :: IO
      INTEGER,DIMENSION(:),INTENT(INOUT) :: DATA
      INTEGER,INTENT(IN) :: LENGTH
      INTEGER,INTENT(INOUT) :: INDEX
      LOGICAL               SWPLUS
      INTEGER               RECORD(128)
!
      DATA       LREC  / 128 /
!
      SWPLUS=.FALSE.
!    GO TO 10
!   ======================================
!      ENTRY DIRAPL(NWF,IO,DATA,LENGTH,INDEX)
!   ======================================
      SWPLUS=.TRUE.
! ----------------------- INITIATE
   10 NREC =INDEX/LREC
      N     =INDEX-LREC*NREC
      NREC =NREC+1
      LMIN  =1
      LMAX  =LREC-N
      IOUSED=IOUSED+LENGTH
      GO TO 100
!     GO TO (40,100),IO
!
!WRITE. (READ FIRST /LAST RECORD)    ------------------- WRITE(IO=1)
!
!  40 IF(N.EQ.0) GO TO 50
!     LMAX  =MIN0(LMAX,LENGTH)
!     GO TO 60
!!
!  50 IF(LMAX.LE.LENGTH) GO TO 70
!     LMAX  =LENGTH
!  60 READ(NWF,REC=NREC,ERR=70) RECORD
!!
!  70 DO 80 L=LMIN,LMAX
!     N     =N+1
!  80 RECORD(N)=DATA(L)
!!---------------------------
!     WRITE(NWF,REC=NREC) RECORD
!!---------------------------
!     IF(LMAX.EQ.LENGTH) GO TO 140
!     N     =0
!     NREC =NREC+1
!     LMIN  =LMAX+1
!     LMAX  =LMAX+LREC
!     GO TO 50
!
!READ-------------------- READ (IO=2)
!
  100 LMAX  =MIN0(LMAX,LENGTH)
!--------------------------
  110 READ(NWF,REC=NREC) RECORD
!--------------------------
      DO 120 L=LMIN,LMAX
      N     =N+1
  120 DATA(L)=RECORD(N)
!
      IF(LMAX.EQ.LENGTH) GO TO 140
      N     =0
      NREC =NREC+1
      LMIN  =LMAX+1
      LMAX  =MIN0(LMAX+LREC,LENGTH)
      GO TO 110
!
  140 IF(SWPLUS) INDEX=INDEX+LENGTH
      RETURN
      END SUBROUTINE DIRACOI

      END MODULE UTILITARIOS_FLOAT

program xsu_read2

    USE UTILITARIOS_FLOAT
    integer,dimension(1) :: NEUTIP,NTHIPO,NOGBIB,MATMAX,NCHANEL,NAXIAL,NDEPEN
    real,dimension(1):: TIMERST
    integer,dimension(1):: ISWREC,ISWCLFD
    integer,dimension(8) :: DATXX
    integer,allocatable :: INDTYP(:),NEUTYP(:),NTHYP(:),NCLVAR(:),NCTYPE(:),IADOTH(:),NFLOTH(:),IADIN(:),NFLIN(:),IADPOT(:),NFLPOT(:),IADBPI(:),NFLBPI(:),IADTFL(:),NFLTFL(:),IADTCL(:),NFLTCL(:),IADDCL(:),NFLDCL(:)
    real,allocatable :: ZONDAT(:,:),OTHER(:),XIN(:),POTDEN(:),BPINT(:),TFUEL(:),TCOOL(:),DCOOL(:)
    real,allocatable :: XSS(:,:),BU(:)!,BETA(:,:)
    real,allocatable :: VM(:)
    integer,dimension(1) :: IADXSO
    integer,dimension(1) :: NFLXSO
    integer,dimension(1) :: SWHEXA,SWFLAT
    integer :: INDICE,I,J,TABLA,ATT,INDX=0
    integer :: karg = 1
	real :: V=1/(7.78681E-07),DT,T
	real :: BETA = 0.0071826,keff = 1.00173
	logical :: EQUILIBRIUM
	
	character(len=32) :: InputXSU,EqBuffer,TimeBuffer
	
	CALL GETARG(1,InputXSU)
    open(unit = 92,access = 'direct', recl = 512 ,file = TRIM(InputXSU))
	CALL GETARG(2,EqBuffer)
	READ(EqBuffer,'(L5)') EQUILIBRIUM
	IF (EQUILIBRIUM) THEN
		WRITE(*,*) EQUILIBRIUM
	ENDIF
	CALL GETARG(3,TimeBuffer)
	READ(TimeBuffer,'(E12.5)') DT
	
	CALL GETARG(4,TimeBuffer)
	READ(TimeBuffer,'(E12.5)') T
	
    INDICE = 0
    CALL DIRAPL(92,2,NEUTIP ,1,INDICE)
    CALL DIRAPL(92,2,NTHIPO ,1,INDICE)
    CALL DIRAPL(92,2,NOGBIB ,1,INDICE)
    CALL DIRAPL(92,2,MATMAX ,1,INDICE)
    CALL DIRAPL(92,2,NCHANEL,1,INDICE)
    CALL DIRAPL(92,2,NAXIAL ,1,INDICE)
    CALL DIRAPL(92,2,NDEPEN ,1,INDICE)
    CALL DIRAPL(92,2,TIMERST,1,INDICE)
    CALL DIRAPL(92,2,ISWREC ,1,INDICE)
    CALL DIRAPL(92,2,ISWCLFD,1,INDICE)
    CALL DIRAPL(92,2,SWHEXA ,1,INDICE)
    CALL DIRAPL(92,2,SWFLAT ,1,INDICE)
    CALL DIRAPL(92,2,DATXX  ,8,INDICE)

    ALLOCATE ( INDTYP(NEUTIP(1)))
    ALLOCATE ( NEUTYP(MATMAX(1)))
    ALLOCATE ( NTHYP(MATMAX(1)))
	ALLOCATE ( NCTYPE(MATMAX(1)))
	ALLOCATE ( NCLVAR(NDEPEN(1)))
	ALLOCATE ( ZONDAT(MATMAX(1),(NDEPEN(1))))
	ALLOCATE (VM(1))
	
    CALL DIRAPL(92,2,INDTYP,NEUTIP(1),INDICE)
    CALL DIRAPL(92,2,NEUTYP,MATMAX(1),INDICE)
    CALL DIRAPL(92,2,NTHYP,MATMAX(1),INDICE)
	CALL DIRAPL(92,2,NCTYPE,MATMAX(1),INDICE)
	WRITE(*,*) 'NEUTYP',NEUTYP
	DO I= 1,(NDEPEN(1))
        CALL DIRAPL(92,2,ZONDAT(:,I),MATMAX(1),INDICE)
    ENDDO
	
	ALLOCATE ( OTHER(100) )
	ALLOCATE ( IADOTH(100) )
	ALLOCATE ( NFLOTH(100) )
!	ATRIBUTO 1	
	ALLOCATE ( XIN(100) )
	ALLOCATE ( IADIN(100) )
	ALLOCATE ( NFLIN(100) )
	
	ALLOCATE ( POTDEN(100) )
	ALLOCATE ( IADPOT(100) )
	ALLOCATE ( NFLPOT(100) )
	
	ALLOCATE ( BPINT(100) )
	ALLOCATE ( IADBPI(100) )
	ALLOCATE ( NFLBPI(100) )
!	ATRIBUTO 2	
	ALLOCATE ( TFUEL(100) )
	ALLOCATE ( IADTFL(100) )
	ALLOCATE ( NFLTFL(100) )
!	ATRIBUTO 3	
	ALLOCATE ( TCOOL(100) )
	ALLOCATE ( IADTCL(100) )
	ALLOCATE ( NFLTCL(100) )
!	ATRIBUTO 4	
	ALLOCATE ( DCOOL(100) )
	ALLOCATE ( IADDCL(100) )
	ALLOCATE ( NFLDCL(100) )
	
	DO TABLA = 1,NEUTIP(1)
		WRITE(*,*) TABLA,INDTYP(TABLA)
		INDICE = INDTYP(TABLA)
		
		CALL DIRAPL(92,2,VM,1,INDICE)
		CALL DIRAPL(92,2,NCLVAR,NDEPEN(1),INDICE)
		CALL DIRAPL(92,2,IADXSO,1,INDICE)
		CALL DIRAPL(92,2,NFLXSO,1,INDICE)
			
		CALL DIRAPL(92,2,OTHER,100,INDICE)
		CALL DIRAPL(92,2,IADOTH,100,INDICE)
		CALL DIRAPL(92,2,NFLOTH,100,INDICE)
		
		CALL DIRAPL(92,2,XIN,100,INDICE)
		CALL DIRAPL(92,2,IADIN,100,INDICE)
		CALL DIRAPL(92,2,NFLIN,100,INDICE)
		
		CALL DIRAPL(92,2,POTDEN,100,INDICE)
		CALL DIRAPL(92,2,IADPOT,100,INDICE)
		CALL DIRAPL(92,2,NFLPOT,100,INDICE)
		
		CALL DIRAPL(92,2,BPINT,100,INDICE)
		CALL DIRAPL(92,2,IADBPI,100,INDICE)
		CALL DIRAPL(92,2,NFLBPI,100,INDICE)
		
		CALL DIRAPL(92,2,TFUEL,100,INDICE)
		CALL DIRAPL(92,2,IADTFL,100,INDICE)
		CALL DIRAPL(92,2,NFLTFL,100,INDICE)
		
		CALL DIRAPL(92,2,TCOOL,100,INDICE)
		CALL DIRAPL(92,2,IADTCL,100,INDICE)
		CALL DIRAPL(92,2,NFLTCL,100,INDICE)
		
		CALL DIRAPL(92,2,DCOOL,100,INDICE)
		CALL DIRAPL(92,2,IADDCL,100,INDICE)
		CALL DIRAPL(92,2,NFLDCL,100,INDICE)
		
		IF (.NOT. ALLOCATED(XSS)) THEN
			ALLOCATE(XSS (NFLXSO(1),(6+NOGBIB(1))*NOGBIB(1)))
			ALLOCATE(BU  (NFLXSO(1)))
		ENDIF
		
		WRITE(*,*) 'IADXSO',IADXSO
		WRITE(*,*) 'NFLXSO',NFLXSO
		WRITE(*,*) NCLVAR(1)
		INDICE = IADXSO(1)
		
		CALL DIRAPL(92,2,BU,NFLXSO(NCLVAR(1)),INDICE)
		WRITE(*,*) BU
		DO IQ = 1,NFLIN(NCLVAR(1))
			INDX = INDICE
			CALL DIRAPL(92,2,XSS(IQ,:),(6+NOGBIB(1))*NOGBIB(1),INDICE)
			WRITE(*,*) XSS(IQ,:)
			IF (.NOT. EQUILIBRIUM) THEN
				WRITE(*,*) 'NOT IN EQUILIBRIUM'
				XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1)) = (1-BETA)*XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1))/keff*1.0015!si se le coloca el pulso, esto va como piña
!				XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1)) = (1-BETA)*XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1))
			ELSE
				WRITE(*,*) 'EQUILIBRIUM'
				XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1)) = (1-BETA)*XSS(IQ,NOGBIB(1)+1:2*NOGBIB(1))/keff
			ENDIF
			IF (DT .NE. 0.00) THEN
				WRITE(*,*) 'DERIVATIVE'
				XSS(IQ,1:NOGBIB(1)) = XSS(IQ,1:NOGBIB(1)) + 1/(V*DT)
			ELSE
				WRITE(*,*) 'NO DERIVATIVE'
			ENDIF
			WRITE(*,*) XSS(IQ,:)
			CALL DIRAPL(92,1,XSS(IQ,:),(6+NOGBIB(1))*NOGBIB(1),INDX)
		ENDDO
	ENDDO
	
    DEALLOCATE (INDTYP)
    DEALLOCATE (NEUTYP)
    DEALLOCATE (NTHYP)
    DEALLOCATE (NCLVAR)
end program xsu_read2

