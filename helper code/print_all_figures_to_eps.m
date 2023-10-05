function [] = print_all_figures_to_eps()
% function [] = print_all_figures_to_eps()
% 
% prints all open figures to .eps files fig01.eps, fig02.eps, etc. in the
% current directory. be careful not to overwrite!
% 
% J. Lucas Mckay, Ph.D. 2011

%for i = 1:max(get(0,'children'))
%	eval(['print -depsc2 -tiff -r300 -f' num2str(i) ' fig' num2str(i,'%02d') '.eps'])
%end

% edited to work in 2015b
% Aiden Payne 2015

for i = 1:length(get(0,'children'))
	%eval(['print -depsc2 -tiff -r300 -f' num2str(i) ' fig' num2str(i,'%02d') '.eps'])
    %print(['-depsc2 -tiff -r300 -f' num2str(i) ' fig' num2str(i,'%02d') '.eps'])
    print(['fig',num2str(i,'%02d'),'.eps'],'-depsc2','-painters', '-tiff', '-r300')
end
