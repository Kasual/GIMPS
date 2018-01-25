# GIMPS
Great Internet Mersenne Prime Search (GIMPS) Finding World Record Primes Since 1996

Welcome to the Great Internet Mersenne Prime Search!

To use this program you must agree to the terms and conditions,
prize rules, etc. at http://mersenne.org/legal/



__Version 29.4 build 7__ - See the latest [WhatsNew.txt](https://www.mersenne.org/download/whatsnew_294b7.txt) file for a full list of changes.

### Highlights of version 29.4 include
- GIMPS has a new sub-project -- finding (probable) prime Mersenne cofactors.
- Like LL tests, PRP tests now support shift counts to aid in running double-checks.
- PRP tests now support a type of low overhead error checking that almost guarantees correct results even on flaky hardware.
- Because PRP tests are highly reliable, we now offer the option to do PRP tests instead of Lucas-Lehmer primality tests.
- For non-base-2 PRP tests, there is a new option to run each iteration twice and rollback if a mismatch occurs.



## Software Source Code
- __If you use GIMPS source code to find Mersenne primes, you must agree to adhere to the [GIMPS free software license agreement](https://www.mersenne.org/legal/#EULA).__ Other than that restriction, you may use this code as you see fit.
- The source code for the program is highly optimized Intel assembly language. There are many more-readable FFT algorithms available on the web and in textbooks. The program is also completely non-portable. If you are curious anyway, you can
- The GIMPS program is very loosely based on C code written by Richard Crandall. Luke Welsh has started a web page that points to Richard Crandall's program and [other available source code](https://www.mersenne.org/download/freeware.php) that you can use to help search for Mersenne primes.



## Software End User License Agreement ("EULA")
- This EULA applies to __all versions of GIMPS Prime95 and MPrime__ software and source code ("Software").
- __Software is free to download and use indefinitely__ on any computer(s) you own or for which you have permission and authority to install and run Software. Software is not export-restricted.
- __To use the Software you agree to be bound by this EULA and the [Terms and Conditions of Use](https://www.mersenne.org/legal/#TCU).__
- GIMPS reserves the right to change this EULA without notice and with reasonable retroactive effect. Last updated 15 October 2008.
- GIMPS not responsible for any damages or losses arising from use of Software. SOFTWARE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.



## Terms and Conditions of Use ("TCU")
- __GIMPS participation is free of charge and open to the public internationally.__
- _Web Site_. The GIMPS "Web Site" is understood to include all __[Mersenne.org](https://www.mersenne.org/)__ Internet web site domains, web content and services, PrimeNet APIs, data, downloads, etc., regardless of means of access.
- _Non-Participants_. __Viewing the Web Site does not necessarily make you a Participant;__ non-Participants are not bound by this TCU.
- _Participation Constitutes Agreement_. __"Participant" is understood to be__ an individual person, or a single individual person designated as the authorized representative of any group, team, organization or legal entity, who personally, or whose computer(s), __accesses and/or communicates with the Web Site to perform, or cause to perform, mathematical calculations that are, or typically, systematically organized by GIMPS__. This includes, but is not restricted to, use of Prime95, MPrime, GLucas, or similar software, automatically over the Internet or using the Web Site (manual testing forms, reports, APIs, etc.), to get data or work assignments or to upload results or other data. __Participation constitutes agreement with the TCU__ by that individual and any group, team, organization or legal entity the Participant represents.
- _Participant Identifiers_. Participant's user ID, password and email address are the primary unique identifiers used by GIMPS to communicate and confirm Participant's identity. Secondary identifiers may include unique computer identifiers known as a "GUID". The "Anonymous" or "ANONYMOUS" user ID is owned by GIMPS, but may be used by Participants who do not wish to be publically identified.
- _Award Rules_. Participant agrees to the [Research Discovery Award Rules](https://www.mersenne.org/legal/#rules).
- _Award Refusal_. Participant may at their option decline any award. Research Discovery Award Rules apply even if an award is unclaimed or optionally declined by the Participant.
- _Data Ownership and Privacy_. GIMPS owns all collected data, and may publish or make available certain limited subset(s) of that data for public reference, excluding personally identifiable data according to the [Privacy Policy](https://www.mersenne.org/legal/#privacy). Examples of published data are stats, standings, charts and other derived charts or tables.
- _Disclaimer_. GIMPS is not responsible for any losses due to Web Site errors, electronic transmission errors, omissions or unauthorized disclosures, failure of any software to correctly find and timely report a new prime number, or any other research discovery, or for someone "poaching" or "stealing" your assignment (performing work on it without it being assigned to them by GIMPS) and subsequently making a discovery.
- _Terms and Conditions of Use Changes_. GIMPS reserves the right to change this TCU without notice. Last updated 15 October 2008.
- _Jurisdiction_. Jurisdiction of law shall be the State of California and the United States of America. Participant agrees to settle all disputes through a good faith effort directly with GIMPS officers and directors, or as a last resort, by third-party binding arbitration through a certified arbitrator of GIMPS' choosing.
