import os
import shutil
from pathlib import Path

# 현재 경로 설정 (스크립트를 데이터 상위 디렉토리에서 실행한다고 가정)
base_path = Path(".")

# UL로 시작하는 모든 디렉토리 탐색
for year_dir in base_path.glob("UL*"):
    if not year_dir.is_dir():
        continue
        
    mc_dir = year_dir / "MC"
    
    # MC 폴더가 존재하는지 확인
    if mc_dir.exists() and mc_dir.is_dir():
        print(f"작업 중: {mc_dir}")
        
        # MC 폴더 내부의 모든 항목(폴더/파일)을 순회
        for item in mc_dir.iterdir():
            target_path = year_dir / item.name
            
            # 목적지에 동일한 이름이 있는지 확인 (충돌 방지)
            if not target_path.exists():
                shutil.move(str(item), str(target_path))
                print(f"  이동 완료: {item.name} -> {year_dir}/")
            else:
                print(f"  주의: {target_path} 가 이미 존재하여 이동하지 않았습니다.")
        
        # MC 폴더가 비어있다면 삭제
        try:
            mc_dir.rmdir()
            print(f"삭제 완료: 빈 {mc_dir} 폴더 제거됨")
        except OSError:
            print(f"알림: {mc_dir} 내부에 아직 파일이 남아있어 폴더를 삭제하지 않았습니다.")

print("\n모든 작업이 완료되었습니다.")
