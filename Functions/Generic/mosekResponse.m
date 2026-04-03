function mosekResponse(res)
  % Read was successful
  if strcmp(res.rcodestr, 'MSK_RES_OK')
      % Check if we have non-error termination code or OK
      if isempty(strfind(res.rcodestr, 'MSK_RES_ERR'))
      
          solsta = strcat('MSK_SOL_STA_', res.sol.itr.solsta);

          if strcmp(solsta, 'MSK_SOL_STA_OPTIMAL')
              fprintf('An optimal interior-point solution is located\n');
              
          elseif strcmp(solsta, 'MSK_SOL_STA_DUAL_INFEASIBLE_CER')
              fprintf('Dual infeasibility certificate found.');
              task.putintparam(iparam.infeas_report_auto, onoffkey.on)
          elseif strcmp(solsta, 'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')
              fprintf('Primal infeasibility certificate found.');
              task.putintparam(iparam.infeas_report_auto, onoffkey.on)

          elseif strcmp(solsta, 'MSK_SOL_STA_UNKNOWN') 
              % The solutions status is unknown. The termination code 
              % indicates why the optimizer terminated prematurely. 
              fprintf('The solution status is unknown.\n');
              fprintf('Termination code: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);  
          else
            fprintf('An unexpected solution status is obtained.');
          end
      
      else
          fprintf('Error during optimization: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);  
      end

  else
      fprintf('Could not read input file, error: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg)
  end

end