#!/usr/bin/env python3
"""
Motif模块测试
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_import():
    """测试导入"""
    try:
        from plantm6a.analysis import analyze_motifs
        print("✓ 导入plantm6a.analysis.analyze_motifs成功")
        return True
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        print("  提示: 需要先安装pysam: pip install pysam")
        return False

def test_generate_3mers():
    """测试生成3-mer列表"""
    try:
        from plantm6a.analysis import generate_all_3mers_with_A
    except ImportError:
        print("✗ 无法导入generate_all_3mers_with_A (可能pysam未安装)")
        return False
    
    all_3mers = generate_all_3mers_with_A()
    
    assert len(all_3mers) == 16, f"应该有16个3-mer，实际有{len(all_3mers)}个"
    assert 'AAC' in all_3mers, "AAC应该在列表中"
    assert 'GAT' in all_3mers, "GAT应该在列表中"
    
    print(f"✓ 生成了{len(all_3mers)}个3-mer motif: {', '.join(all_3mers[:8])}...")
    return True

if __name__ == '__main__':
    print("="*60)
    print("Motif模块测试")
    print("="*60)
    
    tests = [
        ("导入测试", test_import),
        ("3-mer生成测试", test_generate_3mers),
    ]
    
    passed = 0
    for name, test_func in tests:
        print(f"\n{name}...")
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"✗ 测试失败: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*60}")
    print(f"总结: {passed}/{len(tests)} 个测试通过")
    print("="*60)
