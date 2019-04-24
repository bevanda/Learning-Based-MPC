function [invSet,iterations, Piter] = mpt_maxCtrlSet(sysStruct,Options)
%MPT_MAXCTRLSET     Computes the maximal  robust control invariant set C_inf
%                   or the maximal robust attractive set K_inf
%
%  [invSet,iterations] = mpt_maxCtrlSet(sysStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the maximal control invariant set C_inf or the maximal attractive set K_inf
% Assume a system x(k+1)=f(x(k),u(k),w(k)) subject to constraints x(k) \in X and u(k) \in U; 
% Then
%       C_inf = {x \in R^n | \exists u(k) \in U, s.t. x(k) \in X, \forall k >=0}
%       
% The N-step attractive set K_inf is defined as with respect to the target set Phi: 
%       K_inf = {x \in R^n | \exists u(k) \in U, s.t. x(k) \in X, x(N) \in Phi, \forall k >=0}
%
% Here, by default, Phi is an invariant set around the origin. By default N->Infty.
% The function is able to deal with additive and polytopic uncertainty in the system 
% matrices. 
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% sysStruct        - System structure in the sysStruct format
%
%	See the MPT Manual for additional details on the structure format or   
%   consult one of the example systems (e.g. Double_Integator) which were  
%	provided with this package.                                            
%
% Options.y0bounds - add constraints on y0? (1 - yes, 0 - no)
% Options.Kinf     - If set to 1, then K_inf will be computed. If set to 0,
%                    C_inf will be computed. (default is 0)
% Options.maxCtr   - Maximum number of iterations (default is 1000)
%                    (corresponds to N-step attractive set)
% Options.verbose  - Optional: level of verbosity
% Options.scaling  - Scaling the set at each iteration with a parameter 
%                    0 < lambda < 1 guarantees finite time convergence to a
%                    robust invariant subset of the maximal control invariant
%                    set. (Default: Options.scaling = 1) 
% Options.Vconverge - A non-zero value will force the algorithm to break if
%                     relative increase of volume of atractive set at the next
%                     iteration compared to volume of the set at the previous
%                     iteration decreases below this value. E.g.
%                     Options.Vconverge=1 will terminate the procedure if
%                     (Vnew-Vold)/Vold*100 < 1.
%                     NOTE! Currently works only for LTI systems!
%                     NOTE! Value of this option has percents as units!
%                     NOTE! Should only be used if you are computing Kinf set!!! 
% Options.set_limit - If the invariant set has a chebychev redius which is
%                     smaller than this value the iteration is aborted. 
%                     (Default is 1e-3) 
% Options.useprojection - if true, uses projections to obtain feasible set. if
%                         false, feasible sets are obtained by solving a
%                         multi-parametric program.
%                         (Default is true)
% Options
%   .Q, .R, .Tset  - additional problem-depended options
%   .probStruct    - the whole problem structure can be passed as well
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%  invSet          - maximal control invariant (or attractive set)
%  iterations      - number of iterations that were required to compute set
%  Piter           - sets obtained at each iteration
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% "Robust Low Complexity Feedback Control of Constrained Systems", P. Grieder and M. Morari;
%  submitted
%
% see also MPT_ONESTEPCTRL

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

error(nargchk(1,2,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end
if nargin<2,
    Options = [];
end
if(~isfield(Options,'maxCtr'))
    Options.maxCtr=1000;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'Kinf'),
    Options.Kinf=0;
end
if ~isfield(Options,'y0bounds'),
    Options.y0bounds=1;
end
if ~isfield(Options, 'Vconverge')
    Options.Vconverge = 0;
end 

if Options.Vconverge > 0 & Options.Kinf==0,
    fprintf('Computing Kinf set because Options.Vconverge>0...\n');
    Options.Kinf = 1;
end

[nx,nu,ny,ndyn,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

if isfield(Options, 'probStruct'),
    % use user-defined problem structure
    probStruct = Options.probStruct;
else
    %SET DUMMY VALUES
    probStruct.norm=2; %set to 2 to obtain more efficient computation
    probStruct.N=1; 
    probStruct.Q=eye(nx);
    probStruct.R=eye(nu);
    
    if isfield(Options, 'Tset'),
        % set user-specified target set, if given
        if isa(Options.Tset, 'polytope'),
            if isfulldim(Options.Tset),
                probStruct.Tset = Options.Tset;
            end
        end
    end
    if isfield(Options, 'Q'),
        % set user-specified Q penalty, if given
        probStruct.Q = Options.Q;
    end
    if isfield(Options, 'R'),
        % set user-specified R penalty, if given        
        probStruct.R = Options.R;
    end
end

verOpt.verbose=0;
[dummy, probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);

%USE ONLY PROJECTION TO COMPUTE SET
%(no multiparametric programming)
probStruct.subopt_lev=2;
probStruct.y0bounds=Options.y0bounds;

Options.feasset=1; 

Piter = polytope;

if(Options.Kinf) & ~iscell(sysStruct.A)
    %compute maximal attractive set
    if(~isfield(probStruct,'Tset') | ~isfulldim(probStruct.Tset))
        [Matrices]=mpt_constructMatrices(sysStruct,probStruct,Options);  
        probStruct.Tset=Matrices.Pinvset;
    end
    if nargout>2,
        [invSet,iterations,Piter] = mpt_maxCtrlSetLTI(sysStruct,probStruct,Options);
    else
        [invSet,iterations] = mpt_maxCtrlSetLTI(sysStruct,probStruct,Options);
    end
elseif ~iscell(sysStruct.A)
    %compute maximal control invariant set
    Options.includeLQRset=0;
    if nargout>2,
        [invSet,iterations,Piter] = mpt_maxCtrlSetLTI(sysStruct,probStruct,Options);
    else
        [invSet,iterations] = mpt_maxCtrlSetLTI(sysStruct,probStruct,Options);
    end
elseif Options.Kinf,
    invSet = mpt_iterativePWA(sysStruct,probStruct,Options);
    iterations = [];
else
    [invSet,iterations] = mpt_maxCtrlSetPWA(sysStruct,Options);
end
return


%=====================================================================
function [Pfinal,loopCtr,Piter] = mpt_maxCtrlSetLTI(sysStruct, probStruct, Options)
% compute Cinf/Kinf sets for LTI systems

if(~isfield(Options,'maxCtr'))
    Options.maxCtr=1000;
end
if ~isfield(Options,'set_limit'),
    Options.set_limit=1e-3;
end
if ~isfield(Options, 'useprojection'),
    Options.useprojection = 1;
end
if ~isfield(Options,'scaling'),
    Options.scaling=1;
elseif(Options.scaling>1)
    error('Scaling parameter ''Options.scaling'' must be smaller than 1')
elseif(Options.scaling<=0)
    error('Scaling parameter ''Options.scaling'' must be larger than 0')
end

%initialize
loopCtr=0;
PfinalOld = probStruct.Tset;
notconverged=1;

if nargout > 2,
    % also return sets at individual iterations
    Piter = PfinalOld;
end

while notconverged 
    
    if Options.Vconverge > 0,
        % compute volume of old set
        VfinalOld = volume(PfinalOld);
    end         
    
    loopCtr=loopCtr+1;
    if Options.verbose>1,
        fprintf('Iteration %d      \r',loopCtr);
    else
        if mod(loopCtr,20)==0 | loopCtr==1,
            if Options.verbose > -1,
                fprintf('Iteration %d       \r',loopCtr);
            end
        end
    end
    
    % construct the problem - one step solution to previously computed target
    % set
    Options.includeLQRset = 0;
    tmpProbStruct = probStruct;
    tmpProbStruct.Tset = PfinalOld;
    tmpProbStruct.Tconstraint = 2;
    tmpProbStruct.subopt_lev = 0;
    tmpProbStruct.N = 1;
    
    % get matrices of the problem
    Matrices = mpt_constructMatrices(sysStruct,tmpProbStruct,Options); 
    
    % exit if problem is not feasible
    if isinf(Matrices.W),
        error('Problem is infeasible.');
    end
    
    if Options.useprojection,
        % compute feasible set via projection
        
        % polytope is already in reduced representation, because
        % mpt_constructMatrices calls reduce(). therefore we just mark the
        % polytope as being in minimal representation, saving some computational
        % time.
        P=polytope([-Matrices.E Matrices.G], Matrices.W, 1, 1);  
        
        try
            % first try projection, it's faster
            Pfinal=projection(P, (1:size(Matrices.E,2)), Options);
        catch
            % in case of troubles switch to mpQP
            [dd1,dd2,dd3,dd4,Pfinal] = mpt_mpqp(Matrices);
        end
        
        if Options.verbose > 1,
            fprintf('i = %d, nc(P) = %d, nc(Pfinal) = %d\n', loopCtr, nconstr(P), nconstr(Pfinal));
        end
        
    else
        % solve multi-parametric program
        tmpOptions = Options;
        tmpOptions.verbose = -1;
        if tmpProbStruct.norm == 2,
            [Pn,Fi,Gi,activeConstraints,Pfinal]=mpt_mpqp(Matrices, tmpOptions);
        else
            [Pn,Fi,Gi,activeConstraints,Pfinal]=mpt_mplp(Matrices, tmpOptions);
        end
        
        if Options.verbose > 1,
            fprintf('i = %d, nc(P) = %d, nc(Pfinal) = %d\n', loopCtr, nconstr(PfinalOld), nconstr(Pfinal));
        end
    end
    
    % algorithm converges when two consequent feasible sets are equal
    if(isfulldim(PfinalOld) & Options.scaling<1)
        %set is being "shrunk" from the outside in; only applies for computation of Cinf
        notconverged=~ge(Pfinal,PfinalOld,Options);       % Pfinal => PfinalOld
    elseif(isfulldim(PfinalOld) & Options.scaling==1)
        %no scaling, hence convergence is reached when sets are eual
        %this option can apply to the computation of Kinf and Cinf
        notconverged=~eq(Pfinal,PfinalOld,Options);       % Pfinal == PfinalOld
    else
        notconverged=1;
    end
    
    if Options.scaling<1,
        PfinalOld = Pfinal*Options.scaling;
    else
        PfinalOld = Pfinal;
    end
    if loopCtr > Options.maxCtr,
        disp('Maximum number of iterations reached without convergence! Increase value of Options.maxCtr');
        break
    end
    
    if Options.Vconverge > 0,
        % check whether relative increase of volume of Pfinal is bigger than
        % some given value. If not, abort the computation.
        Vfinal = volume(Pfinal);
        diffV = 100*(Vfinal - VfinalOld)/VfinalOld;
        if Options.verbose > 1,
            fprintf('Increase of volume of Kinf: absolute: %e, relative: %.2f %%\n', Vfinal-VfinalOld, diffV);
        end
        if diffV < Options.Vconverge,
            if Options.verbose > -1,
                fprintf('Volume difference (%.2f %%) below tolerance (%.2f %%), aborting computation...\n', ...
                    diffV, Options.Vconverge);
            end
            notconverged = 0;
        end
    end 
    
    if nargout > 2,
        Piter = [Piter Pfinal];
    end
    
    [x,R]=chebyball(PfinalOld);
    if(~isfulldim(PfinalOld) | R<=Options.set_limit)
        disp(['The invariant set for this system (if it exists) has a chebychev radius smaller than Options.set_limit (=' num2str(Options.set_limit) ') !!'])
        disp('Aborting iteration...')
        Pfinal = polytope;
        return
    end
end



%=====================================================================
function [invSet,iterations]=mpt_maxCtrlSetPWA(sysStruct,Options)
% compute maximum controllable set for PWA systems

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);
probStruct.norm=2; %set to 2 to obtain more efficient computation
probStruct.N=1; 
probStruct.Q=eye(nx);
probStruct.R=eye(nu);
verOpt.verbose=0;
probStruct=mpt_verifyProbStruct(probStruct,verOpt);

if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end
origProbStruct = probStruct;

if nargin<3,
    Options = [];
end
if ~isfield(Options,'maxiterations'),
    Options.maxiterations=100;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end

if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);
ssInfo.nx = nx;
ssInfo.nu = nu;
ssInfo.ny = ny;
ssInfo.nPWA = nPWA;
ssInfo.nbool = nbool;
ssInfo.ubool = ubool;
ssInfo.intInfo = intInfo;

emptypoly=mptOptions.emptypoly;

starttime=cputime;
convergedN = -1;

CSstorage = cell(1,intInfo.stacks);
maxPfinal = cell(1,intInfo.stacks);
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.Pfinals = emptypoly;
    CSstorage{ctr}.dynamics = [];
    maxPfinal{ctr} = emptypoly;
    maxPfinalStep{1}{ctr} = emptypoly;
end
emptyCS = CSstorage;
Pn = sysStruct.Pbnd;
for ii=1:nPWA,
    dyn_stack = intInfo.dyns_links(ii,2);
    CSstorage{dyn_stack}.Pfinals = [CSstorage{dyn_stack}.Pfinals Pn];
    CSstorage{dyn_stack}.dynamics = [CSstorage{dyn_stack}.dynamics ii];
end
Step{1} = CSstorage;
for ii=1:intInfo.stacks,
    maxPfinalStep{1}{ii} = CSstorage{ii}.Pfinals;
end

Options.noNoiseOnTset=1;
Options.ispwa=1;

targetMatrices = cell(1,nPWA);
for ii=1:nPWA,
    targetMatrices{ii} = {};
end

keepM = [];
Pbnd = sysStruct.Pbnd;
[ObndA, Obndb] = double(Pbnd);

nTargets = 0;
startTsetTime = cputime;

OptionsXU=Options;
OptionsXU.reduce=0;
OptionsXU.constructset=0;
OptionsXU.reduce_output=0;

startTime = cputime;
feasible = 0;

if probStruct.subopt_lev==2,
    Options.iterative=1;
end

maxxHull = emptypoly;
converged = 0;
for horizon = 2:Options.maxiterations+1,
    fprintf('iteration %d        \r', horizon-1);
    expand=zeros(1,nPWA);
    maxhull{horizon} = emptypoly;
    
    for ii=1:intInfo.stacks,
        PFstorage{ii}.stack = {};
    end
    
    
    CSstorage = Step{horizon-1};
    newCSstorage = emptyCS;
    
    for source_dyn = 1:nPWA,
        AllTsets = [];
        Pstorage = emptypoly;
        source_dyn_stack = intInfo.dyns_links(source_dyn,2);

        tmpProbStruct = probStruct;
        tmpProbStruct.N = 1;
        tmpProbStruct.Tset = emptypoly;
        tmpProbStruct.Tconstraint = 0;

        % even if 1/Inf norm is requested, we use 2-norm here
        % to get faster projection. New matrices with 1/Inf norm are
        % constructed afterwards.
        tmpProbStruct.norm=2;

        cmOptions = Options;
        cmOptions.pwa_index = source_dyn;
        cmOptions.includeLQRset = 0;
        cmOptions.verbose = 0;
        cmOptions.noReduce = 0;
        BMatrices=mpt_constructMatrices(sysStruct,tmpProbStruct,cmOptions);
        W = BMatrices.W;
        if isinf(-BMatrices.W),
            % problem is infeasible
            continue
        end

        % first extract targets regions from dynamics source_dyn

        source_dyn_ind = find(CSstorage{source_dyn_stack}.dynamics==source_dyn);
        PF = CSstorage{source_dyn_stack}.Pfinals(source_dyn_ind);
        PF_prev = PF;
        % now add the rest
        other_dyn_ind = setdiff(1:length(CSstorage{source_dyn_stack}.dynamics),source_dyn_ind);
        PF = CSstorage{source_dyn_stack}.Pfinals(other_dyn_ind);
        PF_prev = [PF_prev PF];
        for ii=setdiff(1:intInfo.stacks,source_dyn_stack)
            PF_prev = [PF_prev CSstorage{ii}.Pfinals];
        end
        % PF_prev now contains all target sets to explore in this iteration

        for reg=1:length(PF_prev)
            Tset = PF_prev(reg);
            
            [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, Tset,nx,nu,source_dyn);
            if ~mfeas,
                continue
            end

            feasible = 1;
            W = Matrices.W;
            G = Matrices.G;
            E = Matrices.E;
            shrinks = 1;
            if isfulldim(maxPfinal{source_dyn_stack}),
                %%do the regiondiff in the lifted XU space to find out, if the currently explored transition will expand our maximum feasible set
                PP=regiondiffXU(Pbnd,maxPfinal{source_dyn_stack},G,W,E,OptionsXU);
                if isfulldim(PP),
                    shrinks=0;
                end
            else
                %%this is the first exploration, force exploration
                shrinks=1;
            end
            if shrinks,
                P=polytope([-E G],W,0,2);
                if isfulldim(P),
                    pOptions = Options;
                    pOptions.noReduce = 1;
                    Pfinal=projection(P,1:size(E,2),pOptions);
                else
                    Pfinal = emptypoly;
                end
                if isfulldim(Pfinal),
                    if ~isminrep(Pfinal),
                        Pfinal = reduce(Pfinal);
                    end
                    newCSstorage{source_dyn_stack}.dynamics = [newCSstorage{source_dyn_stack}.dynamics source_dyn];
                    [newCSstorage{source_dyn_stack}.Pfinals,keep] = sub_mpt_expandlist(Pfinal, newCSstorage{source_dyn_stack}.Pfinals);
                    notkept = find(keep==0);
                    if ~isempty(notkept),
                        newCSstorage{source_dyn_stack}.dynamics(notkept) = [];
                    end
                end
                shrink(source_dyn) = 1;
            end % expands 
        end % region
    end %source_dyn
    
    for ii=1:intInfo.stacks,
        [newCSstorage{ii}.Pfinals,keep] = reduceunion(newCSstorage{ii}.Pfinals,Options);
        notkept = find(keep==0);
        if ~isempty(notkept),
            newCSstorage{ii}.dynamics(notkept) = [];
        end

        nTargets = nTargets + length(newCSstorage{ii}.Pfinals);
        for jj=1:length(newCSstorage{ii}.Pfinals),
            maxPfinal{ii} = sub_mpt_expandlist(newCSstorage{ii}.Pfinals(jj),maxPfinal{ii});
        end
        maxPfinal{ii} = reduceunion(maxPfinal{ii},Options);
        maxPfinalStep{horizon}{ii} = maxPfinal{ii};
    end

    Step{horizon} = newCSstorage;
    oldSets = emptypoly;
    for ii=1:length(Step{horizon-1}),
        oldSets = [oldSets Step{horizon-1}{ii}.Pfinals];
    end
    newSets = emptypoly;
    for ii=1:length(Step{horizon}),
        newSets = [newSets Step{horizon}{ii}.Pfinals];
    end
    
    if oldSets==newSets
        converged = 1;
        break
    end
end

if ~converged
    error('mpt_maxCtrlSet: maximum number of iterations reached without convergence!');
end

iterations = horizon-1;
invSet = newSets;
return


% ===========================================================================================
function [Pun,keep] = sub_mpt_expandlist(Pf,Pu,Options),
% given a polytope Pf and a polyarray Pu, using a subset check removes all fields of Pu
% which are covered by Pf

Options.reduce=0;        % for fast subset check
Options.constructset=0;  % for fast subset check
Options.elementwise=1;   % go elementwise, i.e. Pu(1)<=Pf, Pu(2)<=Pf, ... , Pu(n)<=Pf
expands = 1;
if ~isfulldim(Pu(1)),
    % if Pu is empty, returns Pf
    Pun=Pf;
    keep=1;
    return
end

PuExt=[Pu Pf];
lenPuExt = length(PuExt);
keep = [];
if lenPuExt>1,
    keep=(~le(Pu,Pf,Options))'; % returns indices of polytopes in Pu which are not a subset of Pf
    keep=[keep 1];
else
    keep=1;
end

Pun=PuExt(find(keep~=0));    % the output consists of polytopes which are not a subset of Pf
return