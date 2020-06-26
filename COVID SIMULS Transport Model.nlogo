extensions [ rngs profiler r ]

globals [

  box-x
  box-y
  box-z
  box-u
  anxietyFactor
  NumberInfected
  InfectionChange
  TodayInfections
  YesterdayInfections
  five
  fifteen
  twentyfive
  thirtyfive
  fortyfive
  fiftyfive
  sixtyfive
  seventyfive
  eightyfive
  ninetyfive
  InitialReserves
  AverageContacts
  AverageFinancialContacts
  ScalePhase
  Days
  GlobalR
  CaseFatalityRate
  DeathCount
  DailyCases
  Scaled_Population
  ICUBedsRequired
  scaled_Bed_Capacity
  currentInfections
  eliminationDate
  PotentialContacts
  bluecount
  yellowcount
  redcount
  todayInfected
  cumulativeInfected
  scaledPopulation
  MeanR
  EWInfections
  StudentInfections
;;  meanDaysInfected

  ;; log transform illness period variables
  Illness_PeriodVariance
  M
  BetaillnessPd
  S


  ;; log transform incubation period variables
  Incubation_PeriodVariance
  MInc
  BetaIncubationPd
  SInc


  ;; log transform compliance period variables
  Compliance_PeriodVariance
  MComp
  BetaCompliance
  SComp

]


breed [ simuls simul ]
breed [ hubs hub ]
breed [ medresources medresource ] ;; people living in the city
breed [ packages package ]


directed-link-breed [red-links red-link]


simuls-own [
  timenow ;; the number of days since initial infection
  inICU ;; whether the person is in ICU or not
  R ;; the estimated RNaught of individuals
  contacts ;; the number of contacts the person has made in the model
  IncubationPd ;; the incubation perios of the illness ascribed to the person
  RiskofDeath ;; the overall risk of deth for the person if they contract the illness based on their age
  Pace ;; the speed that pthe person moves around the environment
  RequireICU ;; a measure of whether the person needs ICU or not
  ownIllnessPeriod ;; unique illness period associated with the individual
  ownIncubationPeriod ;; unique incubation pd for the person - related to IncubationPd so can probably be cleaned up - IncubationPd is a legacy var as previously all incubation periods were identical
  ownComplianceWithIsolation ;; unique variable associated with compliance to Isocation of cases if infected
  asymptom ;; whether the person is asymptomatic or not
  personalVirulence ;; the infectivity of the person
  tracked ;; whether the person has been tracked by the health system
  Asymptomaticflag ;; indicator identifying Asymptomatic cases
  hunted ;; has the person been traced using the phoneApp
  haveApp ;; for use in deterimining if the person has downloaded the app
  wearsMask ;; for use in determining if the person wears a face mask
  studentFlag ;; identifies if the person is a student or not
  wearingMask ;; identifies if the person is wearing a mask or not
  currentVirulence ;; current virulence of the person on the day of their infection
  Imported ;; identifies imported cases
  adultsInHousehold ;; counts how many adults in a household for peole under 70
  propensitytotravel;; likelihood that the person want to travel by air to another destination

  ]

hubs-own [
  settled ;; sets the hub in place
  handrailcleaners ;; does someone clean the handrails?
  temperaturescreening ;;
  physicaldistancing ;;
  transferlikelihood ;;

]

patches-own [
  utilisation ;; indicator of whether any people are located on that patch of the environment or not
  destination ;; indicator of whether this location is a place that people might gather
]

medresources-own [
  capacity ;; bed capacity of hospital system
]

to setup
  profiler:start

  rngs:init
   ;; random-seed  100 ;; for use in setting random nuber generator seeds

  clear-all
  ask patches [ set pcolor black ]
  place-box
  draw-box
  ;; illness period estimation using ln transform
  set Illness_Periodvariance se_Illnesspd
  set BetaIllnessPd  ln ( 1 + ( illness_PeriodVariance / illness_period ^ 2))
  set M ( ln illness_period ) - ( BetaillnessPd / 2)
  set S sqrt BetaIllnessPd

    ;; illness period estimation using ln transform
  set Incubation_Periodvariance se_Incubation
  set BetaIncubationPd  ln ( 1 + ( incubation_PeriodVariance / incubation_period ^ 2))
  set MInc ( ln incubation_period ) - ( BetaincubationPd / 2)
  set SInc sqrt BetaIncubationPd

  ask red-links [ set color red ]
;; sets color of patches to black
;;  ask patches [ if pycor = -40 or pycor = 40  [ set pcolor green ] ]
;;  ask n-of 10 patches with [ pcolor = black ] [ set destination 1 ] ;; a beta function for testing locating many people in one place at a single time

 ;; setting up the hospital
  ask n-of 1 patches [ sprout-medresources 1 ]
  ask medresources [ set color white set shape "Health care" set size 5 set xcor 0 set ycor 10 ]
  calculateScaledBedCapacity
  ask medresources [ ask n-of Scaled_Bed_Capacity patches in-radius 5 [ set pcolor white ] ]


 ;; set up people in the environment and allocates characteristics to them

   create-simuls population
  [ move-to one-of patches with [ pcolor  = black and is-box? box-x box-y box-z box-u = false] set size 1 set shape "person" set color 85 set timenow 0 set IncubationPd int ownIncubationPeriod set InICU 0 set R 0
        set personalVirulence random-normal Global_Transmissability 10 set haveApp random 100
        set wearsMask random 100
        set propensitytotravel random 100
        set ownIllnessPeriod ( exp random-normal M S ) ;; log transform of illness period
        set ownIncubationPeriod ( exp random-normal Minc Sinc ) ;;; log transform of incubation period
        ;;set ownComplianceWithIsolation ( exp random-normal Mcomp SComp )  ;; log transform of compliance with isolation
        rngs:init ;; replacing previous log transform with beta distribution
        let stream_id random-float 999
        let seed random-float 999
        rngs:set-seed stream_id seed
        let dist rngs:rnd-beta  stream_id 28 2
        set ownComplianceWithIsolation dist
        set asymptom random 100
        setASFlag

       ]

  create-hubs 9 [ set shape "airplane" set size 10 set settled false set handrailcleaners random 100 set physicaldistancing random-normal 50 10 set temperaturescreening one-of [ 0 100 ]
    set transferlikelihood (handrailcleaners + physicaldistancing + 100 - temperaturescreening)  ]

  set contact_radius 0 ;; sets contact radius of people
  set days 0 ; used to count days since events - currently redundant
  set Quarantine false
  set eliminationDate 0 ; used to identify the date of elimination where no current, unrecovered cases exist
  set Proportion_People_Avoid PPA ;; used to set the proportion of people who are socially distancing
  set Proportion_Time_Avoid PTA ;; used to set the proportion of time that people who are socially distancing are socially distancing (e.g., 85% of people 85% of the time)
  set spatial_distance false
  set case_isolation false
  placehubs
  linkhubs
  sendtravellerstoNZ

  ;; set up initial infected people

  ask n-of ( Current_Cases ) simuls [  set color red set timenow int ownIncubationPeriod - 1
    set timenow int ownIncubationPeriod - 1 set personalVirulence Global_Transmissability set imported 1 move-to one-of hubs with [ label = "NZ"] ]

  ;; setting households up


;;  ask simuls with [ agerange > 18 and agerange <= 60 ] [ set householdUnit random 600 ] ;; allocates adults to a household unit range
;;  ask simuls with [ agerange > 60 ] [ set householdUnit random 400 + 600 ] ;; allocated older adults to household Units that don't include young children or teenagers
;;  ask simuls with [ agerange > 18 and agerange <= 60 ] [ if count simuls with [ householdUnit = [ householdUnit ] of myself ] > 2 [
;;    set householdUnit random 600 ] ]  ;; allocates up to two adults per household
;;  ask simuls with [ agerange < 19 and studentFlag != 1 ] [ set householdUnit [ householdUnit ] of one-of simuls with [ householdUnit <= 600 and agerange > ([ agerange ] of myself + 20) ] set studentFlag 1  ] ;; Identifies students


  ;; allocates children and teenagers to a household where there are adults at least 20 years older than them and there are not more than 2 adults in the house

 ;; resetHouseholdUnit ;; iterates this process
  set tracking false ;; ensures this is set to false each time the model starts
  set link_switch false ;; ensures this is set to false each timme the model starts
  set schoolspolicy false ;; ensures that the schools settings don't begin before the policy trigger starts
  set maskPolicy false ;; that the mask policy doesn't begin before the policy trigger starts
  set assignAppEss false ;; that the assigning the App to EssentialWorkers doesn't begin before the policy trigger starts
  reset-ticks
end

to placehubs
  ask one-of hubs with [ settled = false ]  [ set xcor -30 set ycor 30 set settled true set heading 0 set label "VIC"]
  ask one-of hubs with [ settled = false ]  [ set xcor -30 set ycor 0 set settled true set heading 0 set label "NSW"]
  ask one-of hubs with [ settled = false ]  [ set xcor -30 set ycor -30 set settled true set heading 0 set Label "SA"]
  ask one-of hubs with [ settled = false ]  [ set xcor 0 set ycor 0 set settled true set heading 0 set label "QLD"]
  ask one-of hubs with [ settled = false ]  [ set xcor 0 set ycor -30 set settled true set heading 0 set label "NZ"]
  ask one-of hubs with [ settled = false ]  [ set xcor 0 set ycor 30 set settled true set heading 0 set label "ACT"]
  ask one-of hubs with [ settled = false ]  [ set xcor 30 set ycor -30 set settled true set heading 0 set label "TAS"]
  ask one-of hubs with [ settled = false ]  [ set xcor 30 set ycor 0 set settled true set heading 0 set label "WA"]
  ask one-of hubs with [ settled = false ]  [ set xcor 30 set ycor 30 set settled true set heading 0 set label "OTHER"]

end

to place-box
    set box-x 20
    set box-y 20
    set box-z -20
    set box-u -20
end

to-report is-box? [ x y z u ] ;; patch reporter
  report
    ( abs pxcor = x and abs pycor <= y) or
    ( abs pycor = y and abs pxcor <= x) or
    ( abs pxcor = u and pycor <= z) or
    ( abs pycor = z and pxcor <= u)

end

to draw-box
  ask patches [
    ifelse pxcor = 20 or pxcor = -20 or pycor = 20 or pycor = -20 [ set pcolor orange ] [ set pcolor black ] ;;
  ]
end

;to resethouseholdUnit ;; allocates children to households
;  if schoolsPolicy = true [
;    ask simuls with [ agerange > 18 and agerange <= 60 ] [ if count simuls with [ householdUnit = [ householdUnit ] of myself ] > 2 and 95 > random 100 [
;    set householdUnit random 600 ] ] ;; allows for upo 5% of houses to be sharehouses / care facilities, etc.
;  ask simuls with [ agerange > 60 ] [ if count simuls with [ householdUnit = [ householdUnit ] of myself ] > 2 and 93 < random 100 [
;    set householdUnit [ householdUnit ] of one-of simuls with [ count other simuls with [ householdUnit = [ householdUnit ] of myself ] = 0  ]]];; allows for older people in group homes to make up to 7% of housing units
;  ]
;end

to sendtravellerstoNZ
  ask simuls [ if propensitytotravel > 90 [ move-to one-of hubs with [ label = "NZ"] fd random 10 ] ]
end

to go ;; these funtions get called each time-step
  ask simuls [ bounce move recover settime death isolation reinfect traceme treat Countcontacts respeed hunt checkMask updatepersonalvirulence avoid ] ;;
  ; *current excluded functions for reducing processing resources**
  ask hubs [ transferviahub ]

  globaltreat
  newCases
  CountInfected
  ForwardTime
  setCaseFatalityRate
  countDailyCases
  calculatePopulationScale
  calculateICUBedsRequired
  calculateScaledBedCapacity
  calculateCurrentInfections
  calculateEliminationDate
  assesslinks
  calculatePotentialContacts
  countRed
  countBlue
  countYellow
  scaledownhatch
  calculateYesterdayInfected
  calculateTodayInfected
  calculateScaledPopulation
  calculateMeanR
  countSchoolInfections
  profilerstop

  finish

  tick

end

to bounce  ;; particle procedure
  ;; if we're not about to hit a wall (yellow patch), or if we're already on a
  ;; wall, we don't need to do any further checks
  ;;if [ pcolor] of patch-ahead 1 = orange [ set pace 0 ]
  ;; get the coordinates of the patch we'll be on if we go forward 1
  let new-px int (xcor + dx)
  let new-py int (ycor + dy)
  ;; if hitting left or right wall, reflect heading around x axis
  if ( new-px = box-x ) or
     ( new-px = box-u  ) and permeability <= random 1000
    [ set heading (- heading) ]
  ;; if hitting top or bottom wall, reflect heading around y axis
  if ( new-py = box-y  ) or
     ( new-py = box-z  ) and permeability <= random 1000
    [ set heading (180 - heading) ]

  if [ pxcor ] of patch-here = max-pxcor [ set heading (0 - heading) fd .5 ]
  if [ pycor ] of patch-here = max-pycor [ set heading (180 - heading) fd .5 ]
  if [ pxcor ] of patch-here = min-pxcor [ set heading (0 - heading) fd .5 ]
  if [ pycor ] of patch-here = min-pycor [ set heading (180 - heading) fd .5 ]

end


to move ;; describes the circumstances under which people can move and infect one another

 if color != red or color != black and spatial_Distance = false [ set heading heading + Contact_Radius + random 45 - random 45 fd pace avoidICUs ] ;; contact radius defines how large the dot of contacts for the person is.

  if any? other simuls-here with [ color = red and asymptomaticFlag = 1 and ( currentVirulence * Asymptomatic_Trans ) > random 100 and wearingMask = 0 ] and color = 85  [
    set color red set timenow 0 traceme ] ;; reduces capacity of asymptomatic people to pass on the virus by 1/3

  if any? other simuls-here with [ color = red and asymptomaticFlag = 0 and currentVirulence > random 100 and wearingMask = 0  ] and color = 85  [
    set color red set timenow 0 traceme ] ;; people who are symptomatic pass on the virus at the rate of their personal virulence, which is drawn from population means

  if any? other simuls-here with [ color = red and asymptomaticFlag = 1 and ( currentVirulence * Asymptomatic_Trans ) > random 100 and wearingMask = 1 ] and color = 85 and random 100 > Mask_Efficacy  [
    set color red set timenow 0 traceme ] ;; accounts for a 56% reduction in transfer through mask wearing

  if any? other simuls-here with [ color = red and asymptomaticFlag = 0 and currentVirulence > random 100 and wearingMask = 1 ] and color = 85 and random 100 > Mask_Efficacy  [
    set color red set timenow 0 traceme ] ;; accounts for a 56% reduction in transfer through mask wearing

  if any? other simuls-here with [ color = 85 ] and color = red and Asymptomaticflag = 1 and ( currentVirulence * Asymptomatic_Trans ) > random 100 and wearingMask = 1 and random 100 > Mask_Efficacy
  [ set R R + 1 set GlobalR GlobalR + 1 ]  ;; asymptomatic and wearing mask
  if any? other simuls-here with [ color = 85 ] and color = red and Asymptomaticflag = 0 and currentVirulence  > random 100 and wearingMask = 1 and random 100 >  Mask_Efficacy
  [ set R R + 1 set GlobalR GlobalR + 1 ] ;; symptomatic and wearing mask

  if any? other simuls-here with [ color = 85 ] and color = red and Asymptomaticflag = 1 and ( currentVirulence * Asymptomatic_Trans ) > random 100 and wearingMask = 0
  [ set R R + 1 set GlobalR GlobalR + 1 ] ;; asymptomatic and not wearing mask
  if any? other simuls-here with [ color = 85 ] and color = red and Asymptomaticflag = 0 and currentVirulence  > random 100 and wearingMask = 0
  [ set R R + 1 set GlobalR GlobalR + 1 ] ;; symptomatic and not wearing mask

    ;; these functions reflect thos above but allow the Reff to be measured over the course of the simulation

  if color = red and Case_Isolation = false and ownCompliancewithIsolation < random 100 [ set heading heading + random 90 - random 90 fd pace ] ;; non-compliant people continue to move around the environment unless they are very sick
  if color = red and Quarantine = false [ avoidICUs ] ;; steers people away from the hospital
  if color = black [ move-to one-of MedResources ] ;; hides deceased simuls from remaining simuls, preventing interaction

end

to avoid ;; these are the circustances under which people will interact
    ifelse
      Spatial_Distance = true and Proportion_People_Avoid > random 100 and Proportion_Time_Avoid > random 100  [
      if any? other simuls-here  [ face one-of neighbors with [ count simuls-here = 0 ] ]]
      ;; so, if the social distancing policies are on and you are distancing at this time and you are not part of an age-isolated group and you are not an essentialworker, then if there is anyone near you, move away if you can.
      ;; else...
       [ set heading heading + contact_Radius fd pace avoidICUs  ] ;; otherwise just move wherever you like

end



to settime
  if color = red [ set timenow timenow + 1  ] ;; asks simuls to start counting the days since they became infected and to also possibly die - dying this way currently not implemented but done at the end of the illness period, instead
end

to recover
  if timenow > Illness_Period and color != black  [
    set color yellow set timenow 0 set inICU 0 set requireICU 0  ] ;; if you are not dead at the end of your illness period, then you become recovered and turn yellow and don;t need hospital resources, anymore
end

to reinfect
  if color = yellow and ReinfectionRate > random 100 [ set color 85 ] ;; if you are recovered but suceptible again, you could become reinfected
end

to avoidICUs
  if [ pcolor ] of patch-here = white and InICU = 0 [ move-to min-one-of patches with [ pcolor = black ]  [ distance myself ]  ] ;; makes sure that simulswho have not been sent to hospital stay outside
end

to GlobalTreat ;; send people to quarantine if they have been identified
  let eligiblesimuls simuls with [ color = red and inICU = 0 and ownIncubationPeriod >= Incubation_Period and asymptom >= AsymptomaticPercentage and tracked = 1 ]
  if (count simuls with [ InICU = 1 ]) < (count patches with [ pcolor = white ]) and Quarantine = true and any? eligiblesimuls ;; only symptomatic cases are identified
    [ ask n-of ( count eligiblesimuls * Track_and_Trace_Efficiency )
      eligiblesimuls [
      move-to one-of patches with [ pcolor = white ] set inICU 1 ]]
end

to treat
     if inICU = 1 and color = red [ move-to one-of patches with [ pcolor = white]  ] ; keeps people withint he bunds of the hospital patches and overrides any other movement so they can;t interact with susceptible people
end

to CalculateDailyGrowth ;; calculated the growth in infectes per day
  set YesterdayInfections TodayInfections
  set TodayInfections ( count simuls with [ color = red  and timenow = 1 ] )
  if YesterdayInfections != 0 [set InfectionChange ( TodayInfections / YesterdayInfections ) ]
end

to countcontacts
  if any? other simuls-here with [ color != black ] [ ;; calculates the number of contacts for simuls
    set contacts ( contacts + count other simuls-here ) ]
end

to death ;; calculates death for individuals and adds them to a total for the population

  if Scalephase = 0 and color = red and timenow = int Illness_Period and RiskofDeath > random-float 1  [ set color black set pace 0 set RequireICU 0 set deathcount deathcount + 1 ]
  if Scalephase = 1 and color = red and timenow = int Illness_Period and RiskofDeath > random-float 1  [ set color black set pace 0 set RequireICU 0 set deathcount deathcount + 10 ]
  if Scalephase = 2 and color = red and timenow = int Illness_Period and RiskofDeath > random-float 1  [ set color black set pace 0 set RequireICU 0 set deathcount deathcount + 100 ]
  if Scalephase = 3 and color = red and timenow = int Illness_Period and RiskofDeath > random-float 1  [ set color black set pace 0 set RequireICU 0 set deathcount deathcount + 1000 ]
  if Scalephase = 4 and color = red and timenow = int Illness_Period and RiskofDeath > random-float 1  [ set color black set pace 0 set RequireICU 0 set deathcount deathcount + 10000 ]
end

to respeed
  if tracked != 1 [ set pace speed ] ;; If people aren't tracked they can move as they wish
end

;to checkutilisation
;  ifelse any? simuls-here [ set utilisation 1 ] [ set utilisation 0 ] ;; records which patches are being occupied by simuls
;end

to scaledown ;; preverses the procedure above after the peak of the epidemic
  if scale = true and count simuls with [ color = red ] <= 25 and yellowcount > redcount and days > 0 and scalephase > 0 [ ask n-of (count simuls with [ color = red ] * .9 ) simuls with [ color = red ]
    [ hatch 10 move-to one-of patches with [ pcolor = black ] ]
  set contact_Radius Contact_radius - (90 / 4) set scalephase scalephase - 1  ]
end

to scaledownhatch ;; removes excess simuls fromt the scaled-down view
  if count simuls > Population [  ask n-of ( count simuls - Population ) simuls with [ color = yellow or color = 85 ] [ die ] ]
end

to forwardTime
  set days days + 1 ;; counts days per tick, likely redundant at present as days are not used for anything right now.
end

to CountInfected ;; global infection count
  set numberinfected cumulativeInfected

end

to setCaseFatalityRate ;; calculates death rate per infected person over the course of the pandemic
  if Deathcount > 0 [ set casefatalityrate  ( Deathcount / numberInfected ) ]

end

to countDailyCases ;; sets the day for reporting new cases at 6 (adjustable) days after initial infection, scales up as the population scales
  let casestoday count simuls with [ color = red and int timenow = Case_Reporting_Delay ]

  if Scalephase = 0 [ set dailyCases casestoday ]
  if Scalephase = 1 [ set dailyCases casestoday * 10 ]
  if Scalephase = 2 [ set dailyCases casestoday * 100 ]
  if Scalephase = 3 [ set dailyCases casestoday * 1000 ]
  if Scalephase = 4 [ set dailyCases casestoday * 10000 ]

end

to calculatePopulationScale ;; population scaling function
  if scalephase = 0 [ set Scaled_Population ( count simuls ) ]
  if scalephase = 1 [ set Scaled_Population ( count simuls ) * 10 ]
  if scalephase = 2 [ set Scaled_Population ( count simuls ) * 100 ]
  if scalephase = 3 [ set Scaled_Population ( count simuls ) * 1000 ]
  if scalephase = 4 [ set Scaled_Population ( count simuls ) * 10000 ]
end

to CalculateICUBedsRequired ;; calculates the number of ICU beds required at any time
  let needsICU count simuls with [ color = red and requireICU = 1 ]

  if scalephase = 0 [ set ICUBedsRequired needsICU  ]
  if scalephase = 1 [ set ICUBedsRequired needsICU * 10 ]
  if scalephase = 2 [ set ICUBedsRequired needsICU * 100]
  if scalephase = 3 [ set ICUBedsRequired needsICU * 1000 ]
  if scalephase = 4 [ set ICUBedsRequired needsICU * 10000 ]

end

to calculateScaledBedCapacity ;; scales the number of patches in the environment that represents Australian bed capacity
   set scaled_Bed_Capacity ( Hospital_Beds_In_Australia / 2500 )
end

to calculateCurrentInfections ;; calculates the number of infected people in the population
   let infectedsimuls count simuls with [ color = red ]

   if Scalephase = 0 [ set currentInfections infectedsimuls ]
   if Scalephase = 1 [ set currentInfections infectedsimuls * 10 ]
   if Scalephase = 2 [ set currentInfections infectedsimuls * 100 ]
   if Scalephase = 3 [ set currentInfections infectedsimuls * 1000 ]
   if Scalephase = 4 [ set currentInfections infectedsimuls * 10000 ]

end

to movepackages
  set heading heading + 5 - 5 fd .5 ;; makes stimulus packages drift in the environment

end

to calculateEliminationDate
  if ticks > 1 and count simuls with [ color = red ] = 0 and eliminationDate = 0 [ set eliminationDate ticks stop ] ;; records the day that no infected people remain in the environment
end


;;;;;;;;;;;;;; *****TRACKING AND TRACING FUNCTIONS*********;;;;;;;;;

to traceme
  if tracked != 1 and tracking = true [  if color = red and track_and_trace_efficiency > random-float 1 [ set tracked 1 ] ] ;; this represents the standard tracking and tracing regime
   if color != red and count my-in-links = 0 [ set hunted 0 set tracked 0 ] ;; this ensures that hunted people are tracked but that tracked people are not necessarily hunted

end


to isolation
  if color = red and ownCompliancewithIsolation > random 100 and tracked = 1 [ ;; tracks people and isolates them even if they are pre incubation period
     set pace 0 ]

  ;; this function should enable the observer to track-down contacts of the infected person if that person is either infected or susceptible.
  ;; it enables the user to see how much difference an effective track and trace system might mack to spread

end

to assesslinks ;; this represents the COVID-Safe or other tracing app function
  if link_switch = true and any? simuls with [ color = red and tracked = 1 and haveApp <= App_Uptake ] [ ask simuls with [ color = red and tracked = 1 and haveApp <= App_Uptake ]
    [ if any? other simuls-here [ create-links-with other simuls-here with [ haveapp <= App_Uptake ] ] ] ;; other person must also have the app installed
    ;; asks tracked simuls who have the app to make links to other simuls who also have the app they are in contact with
  ask simuls with [ haveApp <= App_Uptake ] [ ask my-out-links [ set color blue ] ] ;; Covid-safe app out-links  are set to blue
  ask simuls with [ haveApp > App_Uptake ] [ ask my-in-links [ set color red ] ] ;; in-links red but if there is an out and in-link it will be grey

  ask simuls with [ color != red ] [ ask my-out-links [ die ] ] ;; asks all links coming from the infected agent to die
  ask simuls with [ color = yellow ] [ ask my-in-links [ die ] ] ;; asks all links going to the recovered agent to die
  ]

end

to hunt ;; this specifically uses the app to trace people
  if link_switch = true [
    if Track_and_Trace_Efficiency * TTIncrease > random-float 1 and count my-links > 0 and haveApp <= App_Uptake [ set hunted 1 ]  ;; I need to only activate this if the index case is tracked
  if hunted = 1 [ set tracked 1 ]
  ]  ;;

end



;;;;;;;;;;;;*********END OF TTI FUNCTIONS*******;;;;;;;;;;;;;


to calculatePotentialContacts ;; counts the number of people tracked from infected people
   if Scalephase = 0 [ set PotentialContacts ( count links ) ]
   if Scalephase = 1 [ set PotentialContacts ( count links ) * 10 ]
   if Scalephase = 2 [ set PotentialContacts ( count links ) * 100 ]
   if Scalephase = 3 [ set PotentialContacts ( count links ) * 1000 ]
   if Scalephase = 4 [ set PotentialContacts ( count links ) * 10000 ]
end

to countred ;; as per code
  set redCount count simuls with [ color = red ]
end

to countblue ;; as per code
  set blueCount count simuls with [ color = 85 ]
end

to countyellow ;; as per code
  set yellowcount count simuls with [ color = yellow ]
end

to calculateTodayInfected ;; calculates the number of people infected and recorded today for use in conjunction with yesterday's estimate for calculation of daily growth (see below)
  set todayInfected dailycases
end

to calculateYesterdayInfected ;; calculates the number of people infected and recorded today
  set cumulativeInfected cumulativeInfected + todayInfected
end

to calculateScaledPopulation ;; calculates the scaled population for working with smaller environments
  if scalephase = 0 [ set scaledPopulation Total_Population / 10000 ]
  if scalephase = 1 [ set scaledPopulation Total_Population / 1000 ]
  if scalephase = 2 [ set scaledPopulation Total_Population / 100 ]
  if scalephase = 3 [ set scaledPopulation Total_Population / 10 ]
  if scalephase = 4 [ set scaledPopulation Total_Population ]
end

to calculateMeanr
  ifelse any? simuls with [ color = red and timenow = int ownillnessperiod ] [ set meanR ( mean [ R ] of simuls with [ color = red and timenow = int ownillnessperiod ])] [ set MeanR MeanR ] ;; calculates mean Reff for the population
end

to setASFlag
  if asymptom <= asymptomaticPercentage [ set AsymptomaticFlag 1  ] ;;; records an asymptomatic flag for individual people
end

to countSchoolInfections ;; counts infections among school students
   let studentInfects ( count simuls with [ color = red and StudentFlag = 1 ] )
   if Scalephase = 0 [ set studentInfections studentInfects ]
   if Scalephase = 1 [ set studentInfections studentInfects * 10 ]
   if Scalephase = 2 [ set studentInfections studentInfects * 100 ]
   if Scalephase = 3 [ set studentInfections studentInfects * 1000 ]
   if Scalephase = 4 [ set studentInfections studentInfects * 10000 ]
end

to checkMask ;; identifies people who waear a mask
  if maskPolicy = true [
    ifelse wearsMask <= mask_Wearing [ set wearingMask 1 ] [ set wearingMask 0 ] ]

end

to-report meandaysinfected
   report mean [ timenow ] of simuls with [ color = red ]

end

to updatepersonalvirulence ;; creates a triangular distribution of virulence that peaks at the end of the incubation period
  if color = red and timenow <= ownIncubationPeriod [ set currentVirulence ( personalVirulence * ( timenow / ownIncubationPeriod )) ]
  if color = red and timenow > ownIncubationPeriod [ set currentVirulence ( personalVirulence  * ( ( ownIllnessPeriod - timenow ) / ( ownIllnessPeriod - ownIncubationPeriod ))) ]
end

to profilerstop
  if ticks = 30  [
  profiler:stop          ;; stop profiling
  print profiler:report  ;; view the results
    profiler:reset  ]       ;; clear the data
end


;to visitDestination
;  ask simuls [ ;;; sets up destinations where people might gather and set off superspreader events
;    if remainder ticks Visit_Frequency = 0 and any? patches with [ destination = 1 ] in-radius Visit_Radius [ move-to one-of patches with [ destination = 1 ]
;  ]]
;end
 ;; essential workers do not have the same capacity to reduce contact as non-esssential


to newCases
  if mouse-down?  [ ;; lets loose a set of new infected people into the environment
    create-simuls random 50 [ setxy mouse-xcor mouse-ycor set size 1 set shape "person" set color red set timenow random 5 set R 0
        set personalVirulence random-normal Global_Transmissability 10 set haveApp random 100
        set wearsMask random 100

        set ownIllnessPeriod ( exp random-normal M S ) ;; log transform of illness period
        set ownIncubationPeriod ( exp random-normal Minc Sinc ) ;;; log transform of incubation period
        ;;set ownComplianceWithIsolation ( exp random-normal Mcomp SComp )  ;; log transform of compliance with isolation


        rngs:init ;; replacing previous log transform with beta distribution
        let stream_id random-float 999
        let seed random-float 999
        rngs:set-seed stream_id seed
        let dist rngs:rnd-beta  stream_id 28 2
        set ownComplianceWithIsolation dist

        set asymptom random 100
        setASFlag
  ]]
end

to finish
  if count simuls with [ color = red ] = 0 [ stop ]
end

to linkhubs
  ask hubs [ create-link-with one-of other hubs  ]
  ask links [ set color white ]
end

to transferviahub
  let travelers simuls-here with [ propensitytotravel >= 90 ]
  let destinations link-neighbors
  if any? travelers and transferlikelihood > random 1000 [ ask travelers [ move-to one-of destinations ]]
end

to-report percentageinfections
 report count simuls with [ color != 85 ] / count simuls * 100
end

to-report totalsimuls
  report count simuls
end
@#$#@#$#@
GRAPHICS-WINDOW
825
52
1379
607
-1
-1
6.741
1
10
1
1
1
0
0
0
1
-40
40
-40
40
0
0
1
ticks
30.0

BUTTON
103
63
167
97
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
68
108
132
142
Go
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
650
609
768
643
Trace_Patterns
ask n-of 50 simuls with [ color != black ] [ pen-down ] 
NIL
1
T
OBSERVER
NIL
T
NIL
NIL
1

SWITCH
609
23
809
56
spatial_distance
spatial_distance
1
1
-1000

SLIDER
63
158
203
191
Population
Population
0
9000
4000.0
100
1
NIL
HORIZONTAL

SLIDER
63
193
204
226
Speed
Speed
0
5
0.2
.1
1
NIL
HORIZONTAL

SLIDER
609
316
809
349
Illness_period
Illness_period
0
25
20.8
.1
1
NIL
HORIZONTAL

SWITCH
609
61
807
94
case_isolation
case_isolation
1
1
-1000

BUTTON
138
108
202
142
Go Once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
609
353
809
386
ReInfectionRate
ReInfectionRate
0
100
0.0
1
1
NIL
HORIZONTAL

SWITCH
609
205
809
238
quarantine
quarantine
1
1
-1000

SLIDER
609
170
808
203
Track_and_Trace_Efficiency
Track_and_Trace_Efficiency
0
1
0.25
.05
1
NIL
HORIZONTAL

SLIDER
609
98
807
131
Proportion_People_Avoid
Proportion_People_Avoid
0
100
85.0
.5
1
NIL
HORIZONTAL

SLIDER
609
133
808
166
Proportion_Time_Avoid
Proportion_Time_Avoid
0
100
85.0
.5
1
NIL
HORIZONTAL

SLIDER
1434
156
1583
189
Treatment_Benefit
Treatment_Benefit
0
10
4.0
1
1
NIL
HORIZONTAL

MONITOR
259
632
318
677
R0
mean [ R ] of simuls with [ color = red and timenow = int Illness_Period ]
2
1
11

SWITCH
56
473
200
506
policytriggeron
policytriggeron
1
1
-1000

SLIDER
609
243
808
276
Compliance_with_Isolation
Compliance_with_Isolation
0
100
95.0
5
1
NIL
HORIZONTAL

INPUTBOX
56
330
212
391
current_cases
1.0
1
0
Number

INPUTBOX
56
395
212
456
total_population
2.5E7
1
0
Number

SLIDER
609
393
811
426
Incubation_Period
Incubation_Period
0
10
5.1
.1
1
NIL
HORIZONTAL

PLOT
0
632
384
831
Active (red) and Total (blue) Infections ICU Beds (black)
NIL
NIL
0.0
10.0
0.0
200.0
true
false
"" "\n"
PENS
"Current Cases" 1.0 1 -7858858 true "" "plot currentInfections "
"Total Infected" 1.0 0 -13345367 true "" "plot NumberInfected "
"ICU Beds Required" 1.0 0 -16777216 true "" "plot ICUBedsRequired "

PLOT
245
469
547
624
New Infections Per Day
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" "if Scalephase = 1 [ plot count simuls with [ color = red and int timenow = Case_Reporting_Delay ] * 10 ] \nif ScalePhase = 2 [ plot count simuls with [ color = red and int timenow = Case_Reporting_Delay ] * 100 ] \nif ScalePhase = 3 [ plot count simuls with [ color = red and int timenow = Case_Reporting_Delay ] * 1000 ]\nif ScalePhase = 4 [ plot count simuls with [ color = red and int timenow = Case_Reporting_Delay ] * 10000 ]"
PENS
"New Cases" 1.0 1 -5298144 true "" "plot count simuls with [ color = red and timenow = Case_Reporting_Delay ] "

SLIDER
606
436
805
469
Contact_Radius
Contact_Radius
0
180
0.0
1
1
NIL
HORIZONTAL

INPUTBOX
439
105
518
173
ppa
85.0
1
0
Number

INPUTBOX
525
105
610
174
pta
85.0
1
0
Number

TEXTBOX
256
98
432
184
Manually enter the proportion of people who avoid (PPA) and time avoided (PTA) here when using the policy trigger switch
12
0.0
0

SLIDER
606
559
811
592
Hospital_Beds_in_Australia
Hospital_Beds_in_Australia
0
200000
65000.0
5000
1
NIL
HORIZONTAL

SLIDER
1433
195
1586
228
Bed_Capacity
Bed_Capacity
0
20
4.0
1
1
NIL
HORIZONTAL

MONITOR
423
778
487
827
Links
count links / count simuls with [ color = red ]
0
1
12

SWITCH
426
739
540
772
link_switch
link_switch
1
1
-1000

INPUTBOX
1436
410
1591
470
maxv
1.0
1
0
Number

INPUTBOX
1436
480
1591
540
minv
0.0
1
0
Number

INPUTBOX
1436
545
1591
605
phwarnings
0.8
1
0
Number

INPUTBOX
1440
615
1595
675
saliency_of_experience
1.0
1
0
Number

INPUTBOX
1633
349
1788
409
care_attitude
0.5
1
0
Number

INPUTBOX
1633
415
1788
475
self_capacity
0.8
1
0
Number

MONITOR
1663
156
1777
201
Potential contacts
PotentialContacts
0
1
11

MONITOR
433
635
535
680
NIL
numberInfected
17
1
11

PLOT
1826
129
2161
252
Distribution of Illness pd
NIL
NIL
10.0
40.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [ ownIllnessPeriod ] of simuls "

INPUTBOX
1660
210
1816
271
se_illnesspd
4.0
1
0
Number

INPUTBOX
1660
274
1816
335
se_incubation
2.25
1
0
Number

PLOT
1829
250
2167
372
Dist_Incubation
NIL
NIL
0.0
15.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [ ownIncubationPeriod ] of simuls"

INPUTBOX
1433
345
1588
405
initialassociationstrength
0.0
1
0
Number

SLIDER
609
280
811
313
AsymptomaticPercentage
AsymptomaticPercentage
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
606
516
812
549
Global_Transmissability
Global_Transmissability
0
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
248
345
424
378
Essential_Workers
Essential_Workers
0
100
0.0
1
1
NIL
HORIZONTAL

SLIDER
246
380
421
413
Ess_W_Risk_Reduction
Ess_W_Risk_Reduction
0
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
249
308
425
341
App_Uptake
App_Uptake
0
100
0.0
1
1
NIL
HORIZONTAL

SWITCH
252
61
357
94
tracking
tracking
1
1
-1000

SLIDER
370
193
482
226
Mask_Wearing
Mask_Wearing
0
100
47.0
1
1
NIL
HORIZONTAL

SLIDER
248
192
366
225
Mask_Efficacy
Mask_Efficacy
0
100
56.0
1
1
NIL
HORIZONTAL

SWITCH
252
272
374
305
schoolsPolicy
schoolsPolicy
1
1
-1000

SWITCH
428
388
556
421
AssignAppEss
AssignAppEss
1
1
-1000

SLIDER
428
348
556
381
eWAppUptake
eWAppUptake
0
1
0.0
.01
1
NIL
HORIZONTAL

SLIDER
253
20
416
53
TTIncrease
TTIncrease
0
5
2.0
.01
1
NIL
HORIZONTAL

MONITOR
496
778
612
823
Link Proportion
count links with [ color = blue ] / count links with [ color = red ]
1
1
11

MONITOR
1430
45
1562
90
EW Infection %
EWInfections / 2500
1
1
11

MONITOR
1433
95
1566
140
Student Infections %
studentInfections / 2500
1
1
11

SWITCH
379
272
531
305
SchoolPolicyActive
SchoolPolicyActive
1
1
-1000

SLIDER
429
308
561
341
SchoolReturnDate
SchoolReturnDate
0
100
30.0
1
1
NIL
HORIZONTAL

SWITCH
249
231
359
264
MaskPolicy
MaskPolicy
1
1
-1000

SLIDER
433
25
606
58
ResidualCautionPPA
ResidualCautionPPA
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
435
61
608
94
ResidualCautionPTA
ResidualCautionPTA
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
1636
485
1792
518
Case_Reporting_Delay
Case_Reporting_Delay
0
20
6.0
1
1
NIL
HORIZONTAL

PLOT
1638
525
2004
675
R Distribution
NIL
NIL
0.0
15.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [ R ] of simuls with [ color != 85 ] "

SLIDER
1606
49
1779
82
Visit_Frequency
Visit_Frequency
0
100
7.0
1
1
NIL
HORIZONTAL

SLIDER
1608
89
1781
122
Visit_Radius
Visit_Radius
0
10
4.0
1
1
NIL
HORIZONTAL

SLIDER
608
479
810
512
Asymptomatic_Trans
Asymptomatic_Trans
0
1
0.25
.01
1
NIL
HORIZONTAL

MONITOR
432
686
540
731
NIL
currentinfections
17
1
11

BUTTON
73
283
191
317
UnTrace
ask turtles [ pen-up ]
NIL
1
T
OBSERVER
NIL
U
NIL
NIL
1

SWITCH
70
568
175
601
scale
scale
0
1
-1000

SLIDER
1018
14
1191
47
permeability
permeability
0
1000
0.0
1
1
NIL
HORIZONTAL

MONITOR
1083
629
1167
674
% infections
percentageinfections
2
1
11

MONITOR
956
629
1078
674
Mean infected days
meandaysinfected
1
1
11

MONITOR
1176
628
1264
673
Count People
totalsimuls
17
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

bed
false
15
Polygon -1 true true 45 150 45 150 90 210 240 105 195 75 45 150
Rectangle -1 true true 227 105 239 150
Rectangle -1 true true 90 195 106 250
Rectangle -1 true true 45 150 60 195
Polygon -1 true true 106 211 106 211 232 125 228 108 98 193 102 213

bog roll
true
0
Circle -1 true false 13 13 272
Circle -16777216 false false 75 75 150
Circle -16777216 true false 103 103 95
Circle -16777216 false false 59 59 182
Circle -16777216 false false 44 44 212
Circle -16777216 false false 29 29 242

bog roll2
true
0
Circle -1 true false 74 30 146
Rectangle -1 true false 75 102 220 204
Circle -1 true false 74 121 146
Circle -16777216 true false 125 75 44
Circle -16777216 false false 75 28 144

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

box 2
false
0
Polygon -7500403 true true 150 285 270 225 270 90 150 150
Polygon -13791810 true false 150 150 30 90 150 30 270 90
Polygon -13345367 true false 30 90 30 225 150 285 150 150

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

health care
false
15
Circle -1 true true 2 -2 302
Rectangle -2674135 true false 69 122 236 176
Rectangle -2674135 true false 127 66 181 233

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

worker1
true
15
Circle -16777216 true false 96 96 108
Circle -1 true true 108 108 85
Polygon -16777216 true false 120 180 135 195 121 245 107 246 125 190 125 190
Polygon -16777216 true false 181 182 166 197 180 247 194 248 176 192 176 192

worker2
true
15
Circle -16777216 true false 95 94 110
Circle -1 true true 108 107 85
Polygon -16777216 true false 130 197 148 197 149 258 129 258
Polygon -16777216 true false 155 258 174 258 169 191 152 196

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Australia" repetitions="1000" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>count turtles</metric>
    <metric>ticks</metric>
    <metric>numberInfected</metric>
    <metric>deathcount</metric>
    <metric>casefatalityrate</metric>
    <metric>ICUBedsRequired</metric>
    <metric>DailyCases</metric>
    <metric>CurrentInfections</metric>
    <metric>EliminationDate</metric>
    <metric>MeanR</metric>
    <enumeratedValueSet variable="OS_Import_Switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxv">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RestrictedMovement">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days_of_cash_reserves">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_Time_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pta">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cruise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TimeLockDownOff">
      <value value="132"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Track_and_Trace_Efficiency">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Treatment_Benefit">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FearTrigger">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Diffusion_Adjustment">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_population">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Triggerday">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_off">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_incubation">
      <value value="2.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spatial_distance">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Global_Transmissability">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_People_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freewheel">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="self_capacity">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Compliance_with_Isolation">
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Illness_period">
      <value value="20.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stimulus">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="WFH_Capacity">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Bed_Capacity">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReInfectionRate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ppa">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Age_Isolation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Severity_of_illness">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProductionRate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phwarnings">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AsymptomaticPercentage">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean_Individual_Income">
      <value value="55000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="current_cases">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Available_Resources">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="saliency_of_experience">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_illnesspd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Beds_in_Australia">
      <value value="4200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Media_Exposure">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initialassociationstrength">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Superspreaders">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="care_attitude">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact_Radius">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hospital_Beds_in_Australia">
      <value value="65000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Incubation_Period">
      <value value="5.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="case_isolation">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="policytriggeron">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Required">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Speed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ess_W_Risk_reduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Essential_Workers">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tracking">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPPA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPTA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Asymptomatic_Trans">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Proportion">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Post_Proportion">
      <value value="0.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Wuhan" repetitions="100" runMetricsEveryStep="true">
    <setup>setup
set current_cases current_cases + random-normal 20 10
set AsymptomaticPercentage AsymptomaticPercentage + random 10 - random 10
set PPA random 100
set PTA random 100</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>count turtles</metric>
    <metric>ticks</metric>
    <metric>numberInfected</metric>
    <metric>deathcount</metric>
    <metric>casefatalityrate</metric>
    <metric>ICUBedsRequired</metric>
    <metric>DailyCases</metric>
    <metric>CurrentInfections</metric>
    <metric>EliminationDate</metric>
    <metric>MeanR</metric>
    <metric>StudentInfections</metric>
    <metric>EWInfections</metric>
    <metric>count simuls with [ Asymptomaticflag = 1 ]</metric>
    <metric>PPA</metric>
    <metric>PTA</metric>
    <enumeratedValueSet variable="maxv">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RestrictedMovement">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days_of_cash_reserves">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_Time_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pta">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cruise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TimeLockDownOff">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Track_and_Trace_Efficiency">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="App_Uptake">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Treatment_Benefit">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FearTrigger">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Diffusion_Adjustment">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_population">
      <value value="11000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Triggerday">
      <value value="53"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_off">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_incubation">
      <value value="2.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spatial_distance">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Global_Transmissability">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_People_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freewheel">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="self_capacity">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Compliance_with_Isolation">
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Illness_period">
      <value value="20.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stimulus">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="WFH_Capacity">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Bed_Capacity">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReInfectionRate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ppa">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Age_Isolation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Severity_of_illness">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProductionRate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phwarnings">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AsymptomaticPercentage">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean_Individual_Income">
      <value value="55000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="current_cases">
      <value value="86"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Available_Resources">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="saliency_of_experience">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_illnesspd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Beds_in_Australia">
      <value value="4200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Media_Exposure">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initialassociationstrength">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Superspreaders">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="care_attitude">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact_Radius">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hospital_Beds_in_Australia">
      <value value="65000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Incubation_Period">
      <value value="5.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="case_isolation">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="policytriggeron">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Required">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Speed">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ess_W_Risk_reduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Essential_Workers">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tracking">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="schoolsPolicy">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="eWAppUptake">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AssignAppEss">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TTIncrease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SchoolPolicyActive">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maskPolicy">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Proportion">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="NZ new" repetitions="1000" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>count turtles</metric>
    <metric>ticks</metric>
    <metric>numberInfected</metric>
    <metric>deathcount</metric>
    <metric>casefatalityrate</metric>
    <metric>ICUBedsRequired</metric>
    <metric>DailyCases</metric>
    <metric>CurrentInfections</metric>
    <metric>EliminationDate</metric>
    <metric>MeanR</metric>
    <enumeratedValueSet variable="OS_Import_Switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxv">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RestrictedMovement">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days_of_cash_reserves">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_Time_Avoid">
      <value value="89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pta">
      <value value="89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cruise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TimeLockDownOff">
      <value value="99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Track_and_Trace_Efficiency">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Treatment_Benefit">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FearTrigger">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Diffusion_Adjustment">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_population">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Triggerday">
      <value value="39"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_off">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_incubation">
      <value value="2.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spatial_distance">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Global_Transmissability">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_People_Avoid">
      <value value="89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freewheel">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="self_capacity">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Compliance_with_Isolation">
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Illness_period">
      <value value="20.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stimulus">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="WFH_Capacity">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Bed_Capacity">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReInfectionRate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ppa">
      <value value="89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Age_Isolation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Severity_of_illness">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProductionRate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phwarnings">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AsymptomaticPercentage">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean_Individual_Income">
      <value value="55000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="current_cases">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Available_Resources">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="saliency_of_experience">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_illnesspd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Beds_in_Australia">
      <value value="4200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Media_Exposure">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initialassociationstrength">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Superspreaders">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="care_attitude">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact_Radius">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hospital_Beds_in_Australia">
      <value value="65000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Incubation_Period">
      <value value="5.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="case_isolation">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="policytriggeron">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Required">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Speed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ess_W_Risk_reduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Essential_Workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tracking">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPPA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPTA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Asymptomatic_Trans">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Proportion">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Post_Proportion">
      <value value="0.45"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Australia Asymptomatic" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>count turtles</metric>
    <metric>ticks</metric>
    <metric>numberInfected</metric>
    <metric>deathcount</metric>
    <metric>casefatalityrate</metric>
    <metric>ICUBedsRequired</metric>
    <metric>DailyCases</metric>
    <metric>CurrentInfections</metric>
    <metric>EliminationDate</metric>
    <metric>MeanR</metric>
    <enumeratedValueSet variable="OS_Import_Switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxv">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RestrictedMovement">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days_of_cash_reserves">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_Time_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pta">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cruise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TimeLockDownOff">
      <value value="132"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Track_and_Trace_Efficiency">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Treatment_Benefit">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FearTrigger">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Diffusion_Adjustment">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_population">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Triggerday">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_off">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_incubation">
      <value value="2.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spatial_distance">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Global_Transmissability">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_People_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freewheel">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="self_capacity">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Compliance_with_Isolation">
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Illness_period">
      <value value="20.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stimulus">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="WFH_Capacity">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Bed_Capacity">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReInfectionRate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ppa">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Age_Isolation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Severity_of_illness">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProductionRate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phwarnings">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AsymptomaticPercentage">
      <value value="30"/>
      <value value="40"/>
      <value value="50"/>
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean_Individual_Income">
      <value value="55000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="current_cases">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Available_Resources">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="saliency_of_experience">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_illnesspd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Beds_in_Australia">
      <value value="4200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Media_Exposure">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initialassociationstrength">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Superspreaders">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="care_attitude">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact_Radius">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hospital_Beds_in_Australia">
      <value value="65000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Incubation_Period">
      <value value="5.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="case_isolation">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="policytriggeron">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Required">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Speed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ess_W_Risk_reduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Essential_Workers">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tracking">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPPA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ResidualCautionPTA">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Asymptomatic_Trans">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Proportion">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Post_Proportion">
      <value value="0.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Wuhan new" repetitions="100" runMetricsEveryStep="true">
    <setup>setup
set current_cases current_cases + random-normal 20 10
set AsymptomaticPercentage AsymptomaticPercentage + random 10 - random 10
set PPA random 100
set PTA random 100</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>count turtles</metric>
    <metric>ticks</metric>
    <metric>numberInfected</metric>
    <metric>deathcount</metric>
    <metric>casefatalityrate</metric>
    <metric>ICUBedsRequired</metric>
    <metric>DailyCases</metric>
    <metric>CurrentInfections</metric>
    <metric>EliminationDate</metric>
    <metric>MeanR</metric>
    <metric>StudentInfections</metric>
    <metric>EWInfections</metric>
    <metric>count simuls with [ Asymptomaticflag = 1 ]</metric>
    <metric>PPA</metric>
    <metric>PTA</metric>
    <enumeratedValueSet variable="maxv">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RestrictedMovement">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days_of_cash_reserves">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_Time_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pta">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cruise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TimeLockDownOff">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Track_and_Trace_Efficiency">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="App_Uptake">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Treatment_Benefit">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FearTrigger">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Diffusion_Adjustment">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_population">
      <value value="11000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Triggerday">
      <value value="53"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_off">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_incubation">
      <value value="2.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spatial_distance">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Global_Transmissability">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Proportion_People_Avoid">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freewheel">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="self_capacity">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Compliance_with_Isolation">
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Illness_period">
      <value value="20.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stimulus">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="WFH_Capacity">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Bed_Capacity">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReInfectionRate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ppa">
      <value value="85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Age_Isolation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Severity_of_illness">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProductionRate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phwarnings">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AsymptomaticPercentage">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean_Individual_Income">
      <value value="55000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="current_cases">
      <value value="86"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Available_Resources">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="saliency_of_experience">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="se_illnesspd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Beds_in_Australia">
      <value value="4200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Media_Exposure">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initialassociationstrength">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Superspreaders">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="care_attitude">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact_Radius">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hospital_Beds_in_Australia">
      <value value="65000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Incubation_Period">
      <value value="5.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="case_isolation">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="policytriggeron">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ICU_Required">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Speed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ess_W_Risk_reduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Essential_Workers">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tracking">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="schoolsPolicy">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="link_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="eWAppUptake">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AssignAppEss">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="TTIncrease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SchoolPolicyActive">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maskPolicy">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Proportion">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="OS_Import_Switch">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
