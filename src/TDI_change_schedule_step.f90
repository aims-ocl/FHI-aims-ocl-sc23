
!****s* thermodynamic_integration/TDI_change_schedule_step
!  NAME
!    TDI_change_schedule_step
!  SYNOPSIS
subroutine  TDI_change_schedule_step(StepInSchedule)
!  PURPOSE
!     Switch all settings to the new schedule step
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  use thermodynamic_integration
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
  integer,intent(in) :: StepInSchedule
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

   !local
   character*120 :: info_str

   if(use_thermodynamic_integration) then
      TDI_Segment_lattice_vectors(:,:)     = TDI_lattice_vectors(StepInSchedule,:,:)     
      TDI_Segment_force_constants(:,:,:,:) = TDI_force_constants(StepInSchedule,:,:,:,:) 
      TDI_Segment_energy_offset            = TDI_QHA_E0         (StepInSchedule)
   end if 
   TDI_Segment_lambda_start             = TDI_lambda_start   (StepInSchedule)
   TDI_Segment_lambda_end               = TDI_lambda_end     (StepInSchedule)
   TDI_lambda                           = TDI_lambda_start   (StepInSchedule)
   TDI_dlambda_dt                       = (TDI_Segment_lambda_end - TDI_Segment_lambda_start)/MD_schedule_time(StepInSchedule) 
   

   !FIXME In the first step, we do not want to reset the positions after the first cycle
   if (StepInSchedule .gt. 1) then
     if(use_thermodynamic_integration) then
        TDI_Segment_atoms          (:,:)   = TDI_atoms          (StepInSchedule,:,:)     
     end if
     TDI_time_offset                    = TDI_time_offset+MD_schedule_time(StepInSchedule-1)
   end if
  
   ! Info regarding the thermodynamic integration
   write(info_str,'(2X,A,I3,A)') "Thermodynamic integration: Reached end of step",StepInSchedule-1,& 
                                 " in thermodynamic integration schedule"
   call localorb_info(info_str,use_unit,'(A)')
   write(info_str,'(2X,A,I3)')       "| Starting schedule step number :  ",StepInSchedule
   call localorb_info(info_str,use_unit,'(A)')
   if(use_thermodynamic_integration) then
      write(info_str,'(2X,A,E30.15,A)') "| Energy offset                 :  ",TDI_Segment_energy_offset*Hartree," eV"
      call localorb_info(info_str,use_unit,'(A)')
   end if
   write(info_str,'(2X,A,E30.15)')   "| Starting with the lambda value:  ",TDI_Segment_lambda_start
   call localorb_info(info_str,use_unit,'(A)')
   write(info_str,'(2X,A,E30.15)')   "| Stopping   at the lambda value:  ",TDI_Segment_lambda_end
   call localorb_info(info_str,use_unit,'(A)')
   call localorb_info("",use_unit,'(A)')

end subroutine TDI_change_schedule_step
