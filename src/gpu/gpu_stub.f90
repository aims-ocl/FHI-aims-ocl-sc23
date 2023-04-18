subroutine no_gpu()

   use localorb_io, only: localorb_info
   use mpi_tasks, only: aims_stop_coll
   implicit none
   character*120 :: info_str

   write(info_str,'(1X,A)') "* This is a stub GPU call.  You are here because GPU acceleration"
     call localorb_info(info_str)
   write(info_str,'(1X,A)') "* was enabled, but the FHI-aims version you are using was not compiled"
     call localorb_info(info_str)
   write(info_str,'(1X,A)') "* with GPU support. You should not be here. Check the keywords in"
     call localorb_info(info_str)
   write(info_str,'(1X,A)') "* your control.in file or make sure you've compiled aims with GPU support."
     call localorb_info(info_str)
   write(info_str,'(A)') ""
     call localorb_info(info_str)

   call aims_stop_coll( &
     'GPU acceleration was requested, but GPU support not compiled into this FHI-aims &
     &version', 'no_gpu')

end subroutine no_gpu

subroutine set_gpu ()

   use localorb_io
   implicit none
   character*80 :: info_str

   write(info_str,*) "set_gpu."
   call localorb_info(info_str)
   call no_gpu()

end subroutine set_gpu

subroutine get_num_gpus ()

   implicit none

   ! WPH: Skipping debug message here, since I expect end users to actually see
   !      this message.
   call no_gpu()

end subroutine get_num_gpus

subroutine get_gpu_specs ()

   use localorb_io
   implicit none
   character*80 :: info_str

   write(info_str,*) "This is get_gpu_specs."
   call localorb_info(info_str)
   call no_gpu()

end subroutine get_gpu_specs

subroutine initialize_cuda_and_cublas ()

   use localorb_io
   implicit none
   character*80 :: info_str

   write(info_str,*) "This is initialize_cuda_and_cublas."
   call no_gpu ()

end subroutine initialize_cuda_and_cublas

subroutine finalize_cuda_and_cublas ()

   use localorb_io
   implicit none
   character*80 :: info_str

   write(info_str,*) "This is finalize_cuda_and_cublas."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine finalize_cuda_and_cublas

subroutine evaluate_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine evaluate_wave_gpu

subroutine evaluate_hamiltonian_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_hamiltonian_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine evaluate_hamiltonian_shell_gpu

subroutine mgga_contribution_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_mgga_contribution_and_add_to_hamiltonian_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine mgga_contribution_gpu

subroutine add_zora_matrix_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is add_zora_matrix_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine add_zora_matrix_gpu

subroutine set_hamiltonian_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_hamiltonina_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_hamiltonian_gpu

subroutine get_hamiltonian_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_hamiltonina_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_hamiltonian_gpu

subroutine get_hamiltonian_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_hamiltonian_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_hamiltonian_shell_gpu

subroutine hamiltonian_destroy_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is clean_up_variables_integration_hamiltonian_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine hamiltonian_destroy_gpu

subroutine evaluate_ks_density_densmat_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_ks_density_densmat_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine evaluate_ks_density_densmat_gpu

subroutine eval_density_grad_densmat_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_density_gradient_densmat_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_density_grad_densmat_gpu

subroutine update_batch_matrix_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is update_batch_matrix_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine update_batch_matrix_gpu

subroutine hamiltonian_create_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is hamiltonian_create_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine hamiltonian_create_gpu

subroutine density_destroy_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is density_destroy_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine density_destroy_gpu

subroutine update_delta_rho_ks_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is update_delta_rho_ks_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine update_delta_rho_ks_gpu

subroutine update_grad_delta_rho_ks_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is update_grad_delta_rho_ks_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine update_grad_delta_rho_ks_gpu

subroutine set_delta_rho_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_delta_rho_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_delta_rho_gpu

subroutine set_rho_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_rho_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_rho_gpu

subroutine set_rho_change_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_rho_change_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_rho_change_gpu

subroutine set_partition_tab_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_partition_tab_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_partition_tab_gpu

subroutine set_hartree_partition_tab_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_hartree_partition_tab_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_hartree_partition_tab_gpu

subroutine set_delta_rho_gradient_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_delta_rho_gradient_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_delta_rho_gradient_gpu

subroutine set_rho_gradient_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_rho_gradient_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_rho_gradient_gpu

subroutine set_density_matrix_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_density_matrix_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_density_matrix_gpu

subroutine get_delta_rho_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_delta_rho_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_delta_rho_gpu

subroutine get_delta_rho_gradient_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_delta_rho_gradient_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_delta_rho_gradient_gpu

subroutine calculate_rho_change_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is calculate_rho_change_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine calculate_rho_change_gpu

subroutine get_rho_change_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_rho_change_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_rho_change_gpu

subroutine density_create_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is density_create_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine density_create_gpu

subroutine start_profiler_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is start_profiler_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine start_profiler_gpu

subroutine stop_profiler_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is stop_profiler_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine stop_profiler_gpu

subroutine forces_create_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "forces_create_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine forces_create_gpu

subroutine eval_forces_shell_dpsi_h_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_hamiltonian_shell_dpsi_h_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_forces_shell_dpsi_h_psi_gpu

subroutine eval_gga_forces_dens_mat_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_gga_forces_dens_mat_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_gga_forces_dens_mat_gpu

subroutine update_sum_forces_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is update_sum_forces_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine update_sum_forces_gpu

subroutine as_update_sum_forces_and_stress_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is as_update_sum_forces_and_stress_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine as_update_sum_forces_and_stress_gpu

subroutine forces_destroy_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "forces_destroy_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine forces_destroy_gpu

subroutine get_forces_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_forces_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_forces_shell_gpu

subroutine get_sum_forces_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_sum_forces_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_sum_forces_gpu

subroutine get_as_pulay_stress_local_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_as_pulay_stress_local_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_as_pulay_stress_local_gpu

subroutine eval_forces_shell_psi_dh_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is eval_forces_shell_psi_dh_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_forces_shell_psi_dh_psi_gpu

subroutine as_evaluate_gga_stress_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is as_evaluate_gga_stress_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine as_evaluate_gga_stress_gpu

subroutine evaluate_forces_shell_add_mgga_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is evaluate_forces_shell_add_mgga_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine evaluate_forces_shell_add_mgga_gpu

subroutine eval_as_shell_psi_dkin_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is eval_as_shell_psi_dkin_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_as_shell_psi_dkin_psi_gpu

subroutine transpose_as_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is transpose_as_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine transpose_as_shell_gpu

subroutine eval_AS_shell_dpsi_h_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is eval_AS_shell_dpsi_h_psi_gpu"
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_AS_shell_dpsi_h_psi_gpu

subroutine eval_AS_shell_add_psi_kin_psi_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is eval_AS_shell_add_psi_kin_psi_shell_gpu"
   call localorb_info(info_str)
   call no_gpu ()

end subroutine eval_AS_shell_add_psi_kin_psi_shell_gpu

subroutine set_h_times_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_h_times_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_h_times_psi_gpu

subroutine set_d_h_times_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_d_h_times_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_d_h_times_psi_gpu

subroutine set_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_wave_gpu

subroutine set_d_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_d_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_d_wave_gpu

subroutine set_ins_idx_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_ins_idx_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_ins_idx_gpu

subroutine set_permute_compute_by_atom_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_permute_compute_by_atom_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_permute_compute_by_atom_gpu

subroutine set_partition_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_partition_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_partition_gpu

subroutine set_xc_gradient_deriv_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_xc_gradient_deriv_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_xc_gradient_deriv_gpu

subroutine set_xc_tau_deriv_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_xc_tau_deriv_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_xc_tau_deriv_gpu

subroutine set_gradient_basis_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_gradient_basis_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_gradient_basis_wave_gpu

subroutine set_hessian_basis_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_hessian_basis_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_hessian_basis_wave_gpu

subroutine set_dens_mat_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_dens_mat_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_dens_mat_gpu

subroutine set_sum_forces_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_sum_forces_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_sum_forces_gpu

subroutine set_as_pulay_stress_local_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_pulay_stress_local_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_pulay_stress_local_gpu

subroutine set_as_strain_deriv_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_strain_deriv_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_strain_deriv_wave_gpu

subroutine set_as_jac_pot_kin_times_psi_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_jac_pot_kin_times_psi_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_jac_pot_kin_times_psi_gpu

subroutine set_as_hessian_times_xc_deriv_gga_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_hessian_times_xc_deriv_gga_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_hessian_times_xc_deriv_gga_gpu

subroutine set_as_hessian_times_xc_deriv_mgga_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_hessian_times_xc_deriv_mgga_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_hessian_times_xc_deriv_mgga_gpu

subroutine set_as_strain_deriv_kinetic_wave_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_strain_deriv_kinetic_wave_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_strain_deriv_kinetic_wave_gpu

subroutine set_as_strain_deriv_wave_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is set_as_strain_deriv_wave_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine set_as_strain_deriv_wave_shell_gpu

subroutine get_as_strain_deriv_wave_shell_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is get_as_strain_deriv_wave_shell_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine get_as_strain_deriv_wave_shell_gpu

subroutine update_full_matrix_via_map_gpu ()

   use localorb_io
   implicit none

   character*80 :: info_str

   write(info_str,*) "This is update_full_matrix_via_map_gpu."
   call localorb_info(info_str)
   call no_gpu ()

end subroutine update_full_matrix_via_map_gpu

module gpuMR
end module gpuMR

subroutine gpu_MR_stub_base(name)
  use localorb_io, only: localorb_info
  character(*), intent(in) :: name
  call localorb_info('This is '//name)
  call no_gpu()
end subroutine gpu_MR_stub_base

subroutine mr_initialize_gpu()
  call gpu_MR_stub_base('mr_initialize_gpu')
end subroutine mr_initialize_gpu

subroutine evaluate_mr_batch()
  call gpu_MR_stub_base('evaluate_mr_batch')
end subroutine evaluate_mr_batch

subroutine update_mr_batch_gpu()
  call gpu_MR_stub_base('update_mr_batch_gpu')
end subroutine update_mr_batch_gpu

subroutine update_mr_batch_gpu_full()
  call gpu_MR_stub_base('update_mr_batch_gpu_full')
end subroutine update_mr_batch_gpu_full

subroutine get_matrix_packed_gpu()
  call gpu_MR_stub_base('get_matrix_packed_gpu')
end subroutine get_matrix_packed_gpu

subroutine giao_mult_psi_gpu()
  call gpu_MR_stub_base('giao_mult_psi_gpu')
end subroutine giao_mult_psi_gpu

subroutine evaluate_mr_batch_no_symm()
  call gpu_MR_stub_base('evaluate_mr_batch_no_symm')
end subroutine evaluate_mr_batch_no_symm

subroutine symm_antisymm_gpu()
  call gpu_MR_stub_base('symm_antisymm_gpu')
end subroutine symm_antisymm_gpu

subroutine mr_destroy_gpu()
  call gpu_MR_stub_base('mr_destroy_gpu')
end subroutine mr_destroy_gpu
